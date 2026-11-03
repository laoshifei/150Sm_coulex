#!/usr/bin/env python3
"""Run instrumented copies of the original ROOT efficiency macros on source spectra.

Only input selection, output capture and noninteractive display are adapted.
The peak model, initial values, limits, area formula, sqrt(area) errors,
source intensities, 356/344 normalization and Radware fit are retained.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shlex
import shutil
import subprocess


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def replace_function(text, declaration, replacement):
    start = text.index(declaration + '\n{')
    brace = text.index('{', start)
    depth = 1
    end = brace + 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return text[:start] + replacement + text[end:]


PREAMBLE = r'''
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TLatex.h>
#include <TMath.h>
#include <TSystem.h>
#include <Math/MinimizerOptions.h>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;
using json = nlohmann::json;
json result={{"peaks",json::array()}};
double stable_left_tail(double displacement,double delta,double beta,double amplitude) {
  double u=displacement/delta,z=delta/(2.*beta)+u;
  if(z<0.)return amplitude*exp(displacement/beta)*erfc(z);
  double scaled_erfc;
  if(z<25.)scaled_erfc=exp(z*z)*erfc(z);
  else {
    double q=1./(2.*z*z);
    scaled_erfc=(1.-q+3.*q*q-15.*q*q*q+105.*q*q*q*q)/(sqrt(TMath::Pi())*z);
  }
  return amplitude*exp(-u*u-delta*delta/(4.*beta*beta))*scaled_erfc;
}
json fit_json(TF1* f, const TFitResultPtr& r) {
  json p=json::array(),e=json::array(),cov=json::array(),bounds=json::array();
  for(int i=0;i<f->GetNpar();++i) {
    p.push_back(f->GetParameter(i));e.push_back(f->GetParError(i));
    double lo=0,hi=0;f->GetParLimits(i,lo,hi);bounds.push_back({lo,hi});
    json row=json::array();
    for(int j=0;j<f->GetNpar();++j)row.push_back(r.Get()?r->CovMatrix(i,j):0.);
    cov.push_back(row);
  }
  return {{"parameters",p},{"errors",e},{"covariance",cov},{"bounds",bounds},
    {"status",int(r)},{"valid",r.Get()?r->IsValid():false},
    {"covariance_status",r.Get()?r->CovMatrixStatus():-1},
    {"edm",r.Get()?r->Edm():0.},{"chi2",f->GetChisquare()},
    {"ndof",f->GetNDF()},{"fit_points",f->GetNumberFitPoints()}};
}
'''

CAPTURE_PEAK = r'''
  json peak=fit_json(gpeakexregion_tf,peak_result);
  peak["nominal_energy_keV"]=energy;
  peak["area"]=sum;peak["legacy_area_error"]=sqrt(sum);
  peak["range"]={l_bg,l_peak,u_peak,u_bg,ll_ex,lr_ex,hl_ex,hr_ex};
  json bins=json::array();
  for(int ib=1;ib<=h->GetNbinsX();++ib) {
    double x=h->GetBinCenter(ib);
    if(x<l_bg||x>u_bg)continue;
    bool excluded=(x>ll_ex&&x<lr_ex)||(x>hl_ex&&x<hr_ex);
    double peak_value=gstep_tf->Eval(x),background=bgstep_tf->Eval(x);
    bins.push_back({{"energy_keV",x},{"observed",h->GetBinContent(ib)},
      {"peak",peak_value},{"background",background},{"total",peak_value+background},
      {"excluded",excluded}});
  }
  peak["bins"]=bins;result["peaks"].push_back(peak);
'''

RETRY_PEAK = r'''
  json peak_attempts=json::array({fit_json(gpeakexregion_tf,peak_result)});
  if(int(peak_result)!=0) {
    TFitResultPtr best_result=peak_result;
    vector<double> best_par(12),best_err(12);
    for(int ip=0;ip<12;++ip){best_par[ip]=gpeakexregion_tf->GetParameter(ip);best_err[ip]=gpeakexregion_tf->GetParError(ip);}
    double best_chi=gpeakexregion_tf->GetChisquare();int best_ndof=gpeakexregion_tf->GetNDF();
    for(const auto& seed:vector<pair<double,double>>{{.05,1.},{.2,3.}}) {
      gpeakexregion_tf->SetParameters(par);
      gpeakexregion_tf->SetParameter(3,seed.first);gpeakexregion_tf->SetParameter(4,seed.second);
      TFitResultPtr attempt=h->Fit("gpeakexregion_tf","RQS");
      peak_attempts.push_back(fit_json(gpeakexregion_tf,attempt));
      double chi=gpeakexregion_tf->GetChisquare();
      if(isfinite(chi)&&(chi<best_chi-1e-6 || (abs(chi-best_chi)<=1e-6&&int(attempt)==0))) {
        best_chi=chi;best_result=attempt;best_ndof=gpeakexregion_tf->GetNDF();
        for(int ip=0;ip<12;++ip){best_par[ip]=gpeakexregion_tf->GetParameter(ip);best_err[ip]=gpeakexregion_tf->GetParError(ip);}
      }
    }
    gpeakexregion_tf->SetParameters(best_par.data());
    for(int ip=0;ip<12;++ip)gpeakexregion_tf->SetParError(ip,best_err[ip]);
    gpeakexregion_tf->SetChisquare(best_chi);gpeakexregion_tf->SetNDF(best_ndof);
    peak_result=best_result;
  }
'''


def generate(source, ring, stable=False, retry=False):
    text = Path(source).read_text()
    opening = ('TFile *fi = TFile::Open("../rootfile/data0722_build_200ns_jump_300ns.root");\n'
               'TTree *tr = (TTree*)fi->Get("tr_event");\n'
               'TH1D *h = new TH1D("h", "", 4096, 0, 4096);')
    if text.count(opening) != 1:
        raise ValueError('Unexpected original macro input block')
    text = text.replace(opening, 'TH1D *h = nullptr;')
    declaration = 'void get_ring_ge_h(int r)' if ring else 'void get_all_ge_h()'
    text = replace_function(text, declaration, declaration + '\n{ }')
    text = text.replace('h->Fit("gpeakexregion_tf", "RQ");',
                        'TFitResultPtr peak_result=h->Fit("gpeakexregion_tf", "RQS");' +
                        (RETRY_PEAK if retry else ''))
    if stable:
        original = ('Double_t left_skew = gamma*1/2.*a*exp((v[0]-energy)/beta);\n'
                    '  left_skew *= erfc(delta/2./beta + (v[0]-energy)/delta);')
        if text.count(original) != 2:
            raise ValueError('Unexpected historical tail expression')
        text = text.replace(original, 'Double_t left_skew = stable_left_tail(v[0]-energy,delta,beta,gamma*a/2.);')
    text = text.replace('  delete gpeakexregion_tf;', CAPTURE_PEAK + '\n  delete gpeakexregion_tf;')
    if retry:
        text = text.replace('peak["bins"]=bins;', 'peak["attempts"]=peak_attempts;peak["bins"]=bins;')
    text = text.replace('gr->Fit("tf");', 'TFitResultPtr efficiency_result=gr->Fit("tf", "S");')
    text = text.replace('  return tf->GetParameters();', r'''
  result["efficiency"]=fit_json(tf,efficiency_result);
  result["efficiency"]["legacy_residual_rms_percent"]=100.*sqrt(err_all);
  result["efficiency"]["maximum_50_2048"]=max;
  result["efficiency"]["points"]=json::array();
  for(const auto& point:m)result["efficiency"]["points"].push_back(
    {{"energy_keV",point.first},{"value",point.second},{"error",m_err.at(point.first)}});
  return tf->GetParameters();''')
    text = text.replace('  cout << "frac " << frac << endl;',
                        '  result["eu_scale_to_ba"]=frac; cout << "frac " << frac << endl;')
    call = 'eff_cali_ring_ge(std::stoi(argv[3]));' if ring else 'eff_cali_all_ge();'
    main = r'''
int main(int argc,char** argv) {
  if(argc!=4)return 2;
  gROOT->SetBatch(true);
  NUMERICAL_OPTIONS
  unique_ptr<TFile> input(TFile::Open(argv[1],"READ"));
  if(!input||input->IsZombie())return 3;
  h=dynamic_cast<TH1D*>(input->Get(argv[2]));
  if(!h||h->GetNbinsX()!=4096||h->GetBinWidth(1)!=1)return 4;
  result["root_version"]=gROOT->GetVersion();
  result["histogram"]=argv[2];result["histogram_entries"]=h->GetEntries();
  // Original empty-range indexing must never be reached during reconstruction.
  ifstream ranges("fit_range.txt");string line;getline(ranges,line);
  int en=0;double value=0;set<int> energies;
  while(ranges>>en){for(int j=0;j<8;++j)if(!(ranges>>value))return 5;energies.insert(en);}
  for(int e:{40,81,122,245,303,344,356,384,411,444,779,867,964,1112,1408})
    if(!energies.count(e))return 6;
  CALL
  ofstream out("result.json");out<<result.dump(2)<<'\n';
  return out?0:7;
}
'''.replace('CALL', call).replace('NUMERICAL_OPTIONS',
        'ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(30000);' if stable else '')
    return PREAMBLE + '#include <memory>\n#include <set>\n' + text + main


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spectra', type=Path, required=True)
    p.add_argument('--ranges', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--prefix', default='legacy200')
    p.add_argument('--stable-tail', action='store_true',
                   help='Evaluate the same tail stably and allow up to 30000 ROOT function calls; record as numerical recovery')
    p.add_argument('--retry-status', action='store_true',
                   help='For nonzero peak-fit status only, try two fixed alternate tail starts and retain the lowest chi-square')
    p.add_argument('--groups', nargs='+', default=['all', 'ring3', 'ring4', 'ring5'],
                   choices=['all', 'ring3', 'ring4', 'ring5'])
    args = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    spectra, ranges = args.spectra.resolve(strict=True), args.ranges.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    cflags = shlex.split(subprocess.check_output(['root-config', '--cflags'], text=True))
    libs = shlex.split(subprocess.check_output(['root-config', '--libs'], text=True))
    manifest = {'spectra': str(spectra), 'spectra_sha256': sha(spectra),
                'fit_ranges': str(ranges), 'fit_ranges_sha256': sha(ranges),
                'prefix': args.prefix, 'groups': args.groups, 'macros': {},
                'stable_tail_and_30000_calls': args.stable_tail,
                'two_alternate_starts_only_for_nonzero_status': args.retry_status,
                'exact_historical_inputs': False, 'formal_replacement': False}
    for kind in ['all', 'ring']:
        source = root / 'ge_eff' / f'eff_cali_{kind}_ge.cpp'
        adapted = output / f'adapted_{kind}.cpp'
        adapted.write_text(generate(source, kind == 'ring', args.stable_tail, args.retry_status))
        shutil.copyfile(source, output / f'original_{kind}.cpp')
        command = ['g++', '-O2', '-std=c++17', *cflags, str(adapted), '-o',
                   str(output / f'run_{kind}'), *libs]
        with (output / f'compile_{kind}.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        manifest['macros'][kind] = {'path': str(source), 'sha256': sha(source),
                                  'adapted_sha256': sha(adapted), 'compile_command': command}
    for group in args.groups:
        work = output / group
        (work / 'fig').mkdir(parents=True)
        shutil.copyfile(ranges, work / 'fit_range.txt')
        kind = 'all' if group == 'all' else 'ring'
        command = [str(output / f'run_{kind}'), str(spectra), args.prefix + '_' + group,
                   '0' if group == 'all' else group[-1]]
        with (work / 'root_fit.log').open('w') as log:
            subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True,
                           timeout=180)
        result = json.loads((work / 'result.json').read_text())
        print(group, len(result['peaks']), 'peaks;', 'curve status', result['efficiency']['status'],
              'RMS %', result['efficiency']['legacy_residual_rms_percent'], flush=True)
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()

"""Run the real batch script in temporary fresh checkouts with tiny fake executables."""
from pathlib import Path
import json
import os
import re
import subprocess
import tempfile

script_source = (Path(__file__).resolve().parents[2] / 'do_all_urqmd.sh').read_text()
energies = ['3p0', '27']
cut_suffixes = ['', '_strictqLCMS', '_looseqLCMS']
grid_suffixes = ['', '_coarseFFT', '_fineFFT']
fit_suffixes = ['', '_narrowqFit', '_wideqFit']
fit_tags = ['_Cq_fitresults', '_strictqLCMS_Cq_fitresults', '_looseqLCMS_Cq_fitresults',
            '_Cq_fitresults_narrowqFit', '_Cq_fitresults_wideqFit',
            '_Cq_coarseFFT_fitresults', '_Cq_fineFFT_fitresults']
stub = r'''#!/usr/bin/env python3
from pathlib import Path
import json, os, sys, time
base = Path(os.environ['MOCK_BASE'])
name = Path(sys.argv[0]).name
args = sys.argv[1:]
if name == 'sleep':
    time.sleep(.005)
    sys.exit(0)
with (base/'calls.log').open('a') as log:
    log.write(json.dumps([name,args])+'\n')
failure = os.environ.get('MOCK_FAIL','')
prefix = 'UrQMD_3d_source_0-10cent_all_'
if name == 'make':
    for folder in ['analysed','logfiles','figs/fitting/lcms/FFT','figs/syserr/FFT','FFTmethod/results','levyfit/exe','levyfit/object/deps','levyfit/results']:
        assert (base/folder).is_dir(), 'missing fresh-checkout directory: '+folder
    if 'pairsource_urqmd.exe' in args:
        (base/'mode').write_text(next(a[4:] for a in args if a.startswith('FFT=')))
elif name == 'converter_f19.exe':
    Path('AuAu_'+args[0]+'_tree.root').write_text('particles')
elif name == 'pairsource_urqmd.exe':
    if failure == 'analysis': sys.exit(1)
    assert Path('AuAu_'+args[0]+'_tree.root').is_file()
    suffix = ['', '_strictqLCMS', '_looseqLCMS'][int(args[1])]
    mode = '_FFT' if (base/'mode').read_text() == '1' else ''
    output = Path(prefix+args[0]+suffix+mode+'.root')
    output.write_text('writing')
    time.sleep(.02)
    output.write_text('pairs')
elif name == 'transform_fft.exe':
    if failure == 'transform': sys.exit(1)
    source = Path(args[0]); grid = int(args[args.index('--binning')+1])
    assert source.read_text() == 'pairs', 'transform started before analysis completed'
    if grid: assert '_strictqLCMS' not in source.name and '_looseqLCMS' not in source.name
    Path(args[args.index('--output')+1]).write_text('cq')
elif name == 'CqFit.exe':
    if failure == 'fit': sys.exit(1)
    source = Path(args[0]); assert source.read_text() == 'cq'
    variation = int(args[args.index('--fit-range')+1])
    suffix = ['', '_narrowqFit', '_wideqFit'][variation]
    source.with_name(source.name[:-9]+'_fitresults'+suffix+'_FFT.root').write_text('fit')
    if failure == 'invalid_fit' and source.name == prefix+'3p0_Cq_FFT.root' and variation == 0: sys.exit(2)
elif name == 'calc_and_plot_syserr_FFT.exe':
    assert Path(args[args.index('--input-dir')+1]) == base/'analysed'
    selected = args[args.index('--energies')+1].split(',')
    tags = ['_Cq_fitresults', '_strictqLCMS_Cq_fitresults', '_looseqLCMS_Cq_fitresults',
            '_Cq_fitresults_narrowqFit', '_Cq_fitresults_wideqFit', '_Cq_coarseFFT_fitresults', '_Cq_fineFFT_fitresults']
    for energy in selected:
        for tag in tags:
            assert (base/'analysed'/(prefix+energy+tag+'_FFT.root')).read_text() == 'fit', 'final plotter started too early or variant missing'
    (base/'FFTmethod/results/syserr_results_FFT.root').write_text('final')
    if failure == 'final': sys.exit(1)
    if failure == 'final_crash': sys.exit(139)
    if failure in ['invalid_fit','partial_final']: sys.exit(2)
'''


def stem(energy):
    return 'UrQMD_3d_source_0-10cent_all_' + energy


def setup_pairs(base):
    for energy in energies:
        for suffix in cut_suffixes:
            (base/'drho_analyze_urqmd'/(stem(energy)+suffix+'_FFT.root')).write_text('pairs')


def setup_cq(base):
    (base/'analysed').mkdir(exist_ok=True)
    for energy in energies:
        for suffix in cut_suffixes:
            (base/'analysed'/(stem(energy)+suffix+'_Cq_FFT.root')).write_text('cq')
        for suffix in grid_suffixes[1:]:
            (base/'analysed'/(stem(energy)+'_Cq'+suffix+'_FFT.root')).write_text('cq')


def setup_fits(base):
    (base/'analysed').mkdir(exist_ok=True)
    for energy in energies:
        for tag in fit_tags:
            (base/'analysed'/(stem(energy)+tag+'_FFT.root')).write_text('fit')


def run_case(label, settings=None, setup=None, fail='', expected=0):
    base = Path(tempfile.mkdtemp(prefix='fft_batch_'+label+'_', dir='/tmp'))
    for folder in ['bin','drho_analyze_urqmd','FFTmethod','levyfit/exe']:
        (base/folder).mkdir(parents=True, exist_ok=True)
    for executable in ['bin/make','bin/sleep','bin/root.exe','drho_analyze_urqmd/pairsource_urqmd.exe',
                       'drho_analyze_urqmd/converter_f19.exe','FFTmethod/transform_fft.exe','FFTmethod/CqFit.exe',
                       'FFTmethod/calc_and_plot_syserr_FFT.exe','levyfit/exe/EbE_or_Eavg_1d3d_fit.exe']:
        path = base/executable
        path.write_text(stub)
        path.chmod(0o755)
    assert not (base/'figs').exists()
    (base/'rhofitmax_vs_kt_snn.py').write_text('pass\n')
    config = dict(FFTmethod='true',do_conversion='false',do_analysis='true',do_transform='true',do_fitting='true',
                  FFT_systematics='true',do_final_plots='true',FFT_qlcms_cuts='(0 1 2)',
                  energies_low='("3p0")',energies_high='("27")')
    config.update(settings or {})
    script = script_source
    for key,value in config.items():
        script,count = re.subn('^'+key+'=.*$',key+'='+value,script,flags=re.M)
        assert count == 1, key
    (base/'do_all_urqmd.sh').write_text(script)
    for energy in energies:
        (base/'drho_analyze_urqmd'/('AuAu_'+energy+'_tree.root')).write_text('particles')
    old = base/'drho_analyze_urqmd'/(stem('3p0')+'.root')
    old.write_text('original D(rho)')
    if setup: setup(base)
    env = dict(os.environ,MOCK_BASE=str(base),MOCK_FAIL=fail)
    env['PATH'] = str(base/'bin')+':'+env['PATH']
    result = subprocess.run(['bash',str(base/'do_all_urqmd.sh'),'2'],cwd='/tmp',env=env,capture_output=True,text=True,timeout=60)
    (base/'run.log').write_text(result.stdout+result.stderr)
    assert result.returncode == expected, (label,result.returncode,result.stdout,result.stderr,str(base))
    calls = [json.loads(line) for line in (base/'calls.log').read_text().splitlines()]
    if config['FFTmethod'] == 'true':
        assert old.read_text() == 'original D(rho)'
        assert not any(name in ['root.exe','EbE_or_Eavg_1d3d_fit.exe'] for name,args in calls)
    print('PASS batch',label)
    return base,calls


base,calls = run_case('fresh checkout with spaces',{'do_conversion':'true'})
assert sum(name=='pairsource_urqmd.exe' for name,args in calls) == 6
assert sum(name=='transform_fft.exe' for name,args in calls) == 10
assert sum(name=='CqFit.exe' for name,args in calls) == 14
assert sum(name=='calc_and_plot_syserr_FFT.exe' for name,args in calls) == 1

base,calls = run_case('replace_existing',setup=lambda base:(setup_pairs(base),setup_cq(base),setup_fits(base)))
assert sum(name=='pairsource_urqmd.exe' for name,args in calls) == 6
base,calls = run_case('fit_only',{'do_analysis':'false','do_transform':'false'},setup=setup_cq)
assert not any(name in ['pairsource_urqmd.exe','transform_fft.exe'] for name,args in calls)
base,calls = run_case('final_only',{'do_analysis':'false','do_transform':'false','do_fitting':'false'},setup=setup_fits)
assert not any(name in ['pairsource_urqmd.exe','transform_fft.exe','CqFit.exe'] for name,args in calls)
base,calls = run_case('baseline_only',{'FFT_systematics':'false','FFT_qlcms_cuts':'(0)'})
assert sum(name=='CqFit.exe' for name,args in calls) == 2
assert not any(name=='calc_and_plot_syserr_FFT.exe' for name,args in calls)
base,calls = run_case('transform_only',{'do_analysis':'false','do_fitting':'false','do_final_plots':'false'},setup=setup_pairs)
assert sum(name=='transform_fft.exe' for name,args in calls) == 10
assert not any(name in ['CqFit.exe','calc_and_plot_syserr_FFT.exe'] for name,args in calls)
base,calls = run_case('invalid_fit_continues',fail='invalid_fit',expected=2)
assert sum(name=='CqFit.exe' for name,args in calls) == 14
assert any(name=='calc_and_plot_syserr_FFT.exe' for name,args in calls)
for failure in ['analysis','transform','fit','final','final_crash']:
    base,calls = run_case('failure_'+failure,fail=failure,expected=1)
    if failure not in ['final','final_crash']: assert not any(name=='calc_and_plot_syserr_FFT.exe' for name,args in calls)
base,calls = run_case('incomplete_systematics',fail='partial_final',expected=2)
base,calls = run_case('legacy',{'FFTmethod':'false','do_fitting':'false'})
assert not any(name in ['transform_fft.exe','CqFit.exe','calc_and_plot_syserr_FFT.exe'] for name,args in calls)
assert any(name=='make' and 'FFT=0' in args for name,args in calls)
print('All batch checks passed without real event files or production executables.')

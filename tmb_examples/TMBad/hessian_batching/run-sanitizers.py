"""Build and run standalone batching checks against this checkout."""
import argparse,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--eigen-include',required=True);p.add_argument('--output',required=True);p.add_argument('--compiler',default='clang++');a=p.parse_args()
here=Path(__file__).resolve().parent;out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
cmd=[a.compiler,'-std=c++17','-O1','-g','-fno-omit-frame-pointer','-fsanitize=address,undefined','-fno-sanitize-recover=all','-I'+str(here.parents[2]/'TMB/inst/include'),'-I'+a.eigen_include,str(here/'lifecycle.cpp'),'-o',str(out/'lifecycle')]
with (out/'compile.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
env=os.environ.copy();env.update(ASAN_OPTIONS='detect_leaks=0:halt_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
with (out/'run.log').open('w') as f:subprocess.run([str(out/'lifecycle')],env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
print((out/'run.log').read_text().splitlines()[-1])

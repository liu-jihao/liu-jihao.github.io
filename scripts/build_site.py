#!/usr/bin/env python3
"""Generate the static website and editable CV locally. Never publishes."""
from pathlib import Path
import subprocess,sys,argparse
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
for script in ('build_content.py','build_research_topics.py','build_seminar.py','build_discovery.py'):
    cmd=[sys.executable,str(ROOT/'scripts'/script)]+(['--check'] if args.check else [])
    subprocess.run(cmd,check=True,cwd=ROOT)
subprocess.run([sys.executable,str(ROOT/'scripts/build_chrome.py')]+(['--check'] if args.check else []),check=True,cwd=ROOT)
subprocess.run([sys.executable,str(ROOT/'scripts/check_site.py')],check=True,cwd=ROOT)

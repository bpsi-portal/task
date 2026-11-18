#!/usr/bin/env python3
"""
check_dealloc.py - Fortran ALLOCATE/DEALLOCATE 対応チェックスクリプト

使い方:
  python3 check_dealloc.py [ディレクトリ] [オプション]

指定ディレクトリ配下の Fortran ファイルを走査し、ALLOCATE に対応する
DEALLOCATE が同一ディレクトリ内に存在しない変数を検出します。

オプション:
  --missing-only  DEALLOCATE が見つからないものだけ表示
  --summary       パッケージ毎のサマリのみ表示
  --f90-only      .f90 ファイルのみ対象 (.f を除外)
  -h, --help      このヘルプを表示
"""

import os
import re
import sys
import glob


def extract_alloc_vars(line):
    """ALLOCATE文から変数名を抽出する"""
    stripped = line.strip()
    if stripped.startswith('!'):
        return []
    code_part = line.split('!')[0]
    if 'DEALLOCATE' in code_part.upper():
        return []

    m = re.search(r'ALLOCATE\s*\((.+)', code_part, re.IGNORECASE)
    if not m:
        return []

    content = m.group(1)
    content = re.split(r',\s*STAT\s*=', content, flags=re.IGNORECASE)[0]

    vars_found = []
    depth = 0
    current = ''
    for ch in content:
        if ch == '(':
            if depth == 0:
                var_name = current.strip()
                if var_name:
                    vars_found.append(var_name)
                current = ''
            depth += 1
        elif ch == ')':
            depth -= 1
            if depth < 0:
                break
            current = ''
        elif ch == ',' and depth == 0:
            current = ''
        else:
            if depth == 0:
                current += ch

    result = []
    for v in vars_found:
        v = v.strip()
        if not v or v[0].isdigit():
            continue
        if v.upper() in ('STAT', 'IERR', 'SOURCE', 'MOLD'):
            continue
        result.append(v)
    return result


def find_deallocate(var_name, all_lines_by_file):
    """全ファイルの行から DEALLOCATE を検索"""
    search_patterns = [var_name.upper()]
    if '%' in var_name:
        search_patterns.append(var_name.split('%')[-1].upper())

    for fpath, lines in all_lines_by_file.items():
        for lineno, line in lines:
            code_upper = line.split('!')[0].upper()
            if 'DEALLOCATE' not in code_upper:
                continue
            for pat in search_patterns:
                if re.search(r'\b' + re.escape(pat) + r'\b', code_upper):
                    return (fpath, lineno)
    return None


def scan_package(pkg_dir, extensions):
    """1つのパッケージディレクトリをスキャン"""
    # 全ファイルの行を読み込み
    all_lines = {}  # {relpath: [(lineno, line), ...]}
    for root, dirs, files in os.walk(pkg_dir):
        for fname in sorted(files):
            if not any(fname.endswith(ext) for ext in extensions):
                continue
            fpath = os.path.join(root, fname)
            relpath = os.path.relpath(fpath, pkg_dir)
            try:
                with open(fpath, 'r', errors='replace') as f:
                    lines = [(i, l.rstrip()) for i, l in enumerate(f, 1)]
                all_lines[relpath] = lines
            except Exception:
                continue

    if not all_lines:
        return []

    # 継続行を結合
    all_joined = {}
    for relpath, lines in all_lines.items():
        joined = []
        current = ''
        current_no = 0
        for lineno, line in lines:
            if current:
                cont = line.lstrip().lstrip('&').lstrip()
                if line.rstrip().endswith('&'):
                    # continuation line that itself continues
                    cont = cont.rstrip().rstrip('&')
                    current += ' ' + cont
                else:
                    current += ' ' + cont
                    joined.append((current_no, current))
                    current = ''
            else:
                if line.rstrip().endswith('&'):
                    current = line.rstrip()[:-1]
                    current_no = lineno
                else:
                    joined.append((lineno, line))
        if current:
            joined.append((current_no, current))
        all_joined[relpath] = joined

    # ALLOCATE を検出
    results = []
    for relpath, joined in all_joined.items():
        for lineno, line in joined:
            if 'ALLOCATE' not in line.upper():
                continue
            vars_found = extract_alloc_vars(line)
            for var in vars_found:
                dealloc = find_deallocate(var, all_joined)
                if dealloc:
                    results.append(('OK', var, relpath, lineno, dealloc))
                else:
                    results.append(('MISSING', var, relpath, lineno, None))

    return results


def main():
    args = sys.argv[1:]
    target_dir = '.'
    missing_only = False
    summary_only = False
    f90_only = False

    positional = []
    for arg in args:
        if arg == '--missing-only':
            missing_only = True
        elif arg == '--summary':
            summary_only = True
        elif arg == '--f90-only':
            f90_only = True
        elif arg in ('--help', '-h'):
            print(__doc__)
            sys.exit(0)
        else:
            positional.append(arg)

    if positional:
        target_dir = positional[0]

    extensions = ['.f90'] if f90_only else ['.f90', '.f']

    if not os.path.isdir(target_dir):
        print(f"エラー: '{target_dir}' はディレクトリではありません", file=sys.stderr)
        sys.exit(1)

    # ソースファイルが直接あるか確認
    has_direct = any(
        glob.glob(os.path.join(target_dir, f'*{ext}'))
        for ext in extensions
    )

    # サブディレクトリ一覧
    subdirs = sorted([
        d for d in os.listdir(target_dir)
        if os.path.isdir(os.path.join(target_dir, d))
        and not d.startswith('.')
    ])

    # スキャン対象の決定
    if has_direct:
        # ソースが直接ある = 単一パッケージとして扱う
        scan_targets = [('.', target_dir)]
    elif subdirs:
        # サブディレクトリ = 複数パッケージ
        scan_targets = [(d, os.path.join(target_dir, d)) for d in subdirs]
    else:
        print("スキャン対象のファイルが見つかりません")
        sys.exit(0)

    # 色
    C = sys.stdout.isatty()
    RED = '\033[0;31m' if C else ''
    GREEN = '\033[0;32m' if C else ''
    NC = '\033[0m' if C else ''

    print("=" * 60)
    print(" Fortran ALLOCATE/DEALLOCATE 対応チェック")
    print(f" 対象: {os.path.abspath(target_dir)}")
    print(f" 拡張子: {', '.join(extensions)}")
    print("=" * 60)

    total_alloc = 0
    total_missing = 0
    total_ok = 0
    pkg_summary = []

    for pkg_name, pkg_path in scan_targets:
        results = scan_package(pkg_path, extensions)
        if not results:
            continue

        n_ok = sum(1 for r in results if r[0] == 'OK')
        n_miss = sum(1 for r in results if r[0] == 'MISSING')

        total_alloc += len(results)
        total_ok += n_ok
        total_missing += n_miss
        pkg_summary.append((pkg_name, len(results), n_ok, n_miss))

        if summary_only:
            continue

        if missing_only and n_miss == 0:
            continue

        print(f"\n--- {pkg_name} ({n_miss} 漏れ / {len(results)} 総数) ---")
        for status, var, afile, aline, dinfo in results:
            if status == 'MISSING':
                print(f"  {RED}[漏れ]{NC}  {var:<30s} ({afile}:{aline})")
            elif not missing_only:
                df, dl = dinfo
                print(f"  {GREEN}[OK]{NC}    {var:<30s} ({afile}:{aline} -> {df}:{dl})")

    # サマリ
    print("\n" + "=" * 60)
    if len(pkg_summary) > 1:
        print(" パッケージ別サマリ")
        print("=" * 60)
        print(f"  {'パッケージ':<20s} {'総数':>6s} {'OK':>6s} {'漏れ':>6s}")
        print("  " + "-" * 42)
        for pkg, nt, nok, nm in pkg_summary:
            mark = f"{RED}***{NC}" if nm > 0 else ""
            print(f"  {pkg:<20s} {nt:>6d} {nok:>6d} {nm:>6d} {mark}")
        print()

    print(f" ALLOCATE 総数:     {total_alloc}")
    print(f" {GREEN}DEALLOCATE あり:    {total_ok}{NC}")
    print(f" {RED}DEALLOCATE なし:    {total_missing}{NC}")
    if total_alloc > 0:
        rate = (total_ok * 100) // total_alloc
        print(f" カバー率:          {rate}%")
    print("=" * 60)

    sys.exit(1 if total_missing > 0 else 0)


if __name__ == '__main__':
    main()

import os
import time
import subprocess

test_dir = os.path.dirname(os.path.abspath(__file__))

unit_tests = [f for f in os.listdir(test_dir) if f.startswith('test-') and f.endswith('.py')]
n_tests = len(unit_tests)


if __name__ == "__main__":

    print('\n', f'Running all tests → found {n_tests} unit tests:', '\n')

    for it,test in enumerate(unit_tests[:-1]):

        test_path = os.path.join(test_dir, test)

        print('=' * 40)

        print(f'Start test {it+1}/{n_tests} — "{test}"...', '\n')

        t0 = time.time()
        result = subprocess.run(['py', test_path], capture_output=True, text=True)
        t1 = time.time()
        
        if result.returncode == 0:
            print(f'"{test}" passed. — Execution time: {t1 - t0:.2f} seconds', '\n')
        else:
            print(f'"{test} failed with error:\n{result.stderr}', '\n')

        print('=' * 40, '\n')
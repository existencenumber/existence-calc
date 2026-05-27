name: Run Collapse Monster Tests

on:
  push:
    branches: [ main, master ]
  pull_request:
    branches: [ main, master ]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
    - name: Checkout code
      uses: actions/checkout@v4

    - name: Set up Python
      uses: actions/setup-python@v5
      with:
        python-version: '3.10'

    - name: Install dependencies
      run: |
        pip install -r requirements.txt
        pip install requests

    - name: Run unit tests
      run: |
        python tests/test_core.py

    - name: Run API tests
      run: |
        # 启动服务
        python app.py &
        sleep 5
        # 测试各个端点
        python -c "
import requests, json, sys
tests = [
    ('Sum(n**2,(n,1,oo))', 0.0),
    ('Sum(n,(n,1,oo))', -1/12),
    ('Sum(2**n,(n,0,oo))', -1.0),
    ('Sum(factorial(n),(n,0,oo))', 0.5963473623231941),
    ('Sum((-1)**(n+1)/n,(n,1,oo))', 0.6931471805599453),
    ('Sum(n**n,(n,1,oo))', -1.038),  # 理论值
]
failed = False
for query, expected in tests:
    resp = requests.post('http://127.0.0.1:5000/api/calc', json={'query': query})
    data = resp.json()
    if data['status'] != 'success':
        print(f'FAIL: {query} status={data[\"status\"]}')
        failed = True
        continue
    val = data['value']
    # 允许相对误差 1%
    if not abs(val - expected) <= 0.01 * abs(expected):
        print(f'FAIL: {query} value={val} expected={expected}')
        failed = True
    else:
        print(f'PASS: {query} value={val}')
if failed:
    sys.exit(1)
"

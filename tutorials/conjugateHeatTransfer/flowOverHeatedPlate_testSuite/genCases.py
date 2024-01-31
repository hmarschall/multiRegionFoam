import casefoam

# Specify the base case name and a directory name for the tests
baseCase = 'baseCase'
baseCaseDir = 'tests'

# Specify the names of the test cases
caseStructure = [  # main cases
    ['monolithic', 'partitioned_Aitken'],
    # subcases for each of the main cases
    ['Pr_0.01_Re_500_k_1', 'Pr_0.01_Re_500_k_5', 'Pr_0.01_Re_500_k_20',
     'Pr_0.01_Re_10000_k_1', 'Pr_0.01_Re_10000_k_5', 'Pr_0.01_Re_10000_k_20',
     'Pr_100_Re_500_k_1', 'Pr_100_Re_500_k_5', 'Pr_100_Re_500_k_20']]


# Partitioned coupling with acceleration type and initial relaxation factor
def p_coupled(accType, relaxValue):
    return {
        '0/fluid/orig/partitioned/T': {'boundaryField':
                                       {'interface':
                                        {'accType': accType,
                                         'relax': relaxValue}}},

        '#!bash': 'for dir in tests/partitioned_*/*; do cp "baseCase/AllrunP" "$dir/Allrun"; done'
    }


# Monolithic coupling
m_coupled = {
    '#!bash': 'for dir in tests/monolithic/*; do cp "baseCase/AllrunM" "$dir/Allrun"; done'}


# Parameters
def update_params(Pr, Re, k):
    L = 1.0
    rhof = 1.0
    ks = 100.0
    Uinf = 1.0
    mu = rhof * Uinf * L / int(Re)
    kf = ks / k
    cp = kf * Pr / mu
    return {
        'constant/fluid/transportProperties':
        {
            'mu': 'mu [ 1 -1 -1 0 0 0 0 ] %s' % mu,
            'cp': 'cp [ 0 2 -2 -1 0 0 0 ] %s' % cp,
            'k': 'k [1 1 -3 -1 0 0 0] %s' % kf
        },

        '0/fluid/orig/monolithic/k': {'internalField': 'uniform %s' % kf}
    }


# Specify parameters for each of the cases listed in "caseStructure"
caseData = {
    'monolithic': m_coupled,
    'partitioned_Aitken': p_coupled('aitken', '0.75'),
    'Pr_0.01_Re_500_k_1': update_params(0.01, 500, 1),
    'Pr_0.01_Re_500_k_5': update_params(0.01, 500, 5),
    'Pr_0.01_Re_500_k_20': update_params(0.01, 500, 20),
    'Pr_0.01_Re_10000_k_1': update_params(0.01, 10000, 1),
    'Pr_0.01_Re_10000_k_5': update_params(0.01, 10000, 5),
    'Pr_0.01_Re_10000_k_20': update_params(0.01, 10000, 20),
    'Pr_100_Re_500_k_1': update_params(100, 500, 1),
    'Pr_100_Re_500_k_5': update_params(100, 500, 5),
    'Pr_100_Re_500_k_20': update_params(100, 500, 20),
}

# Generate 'tests' directory
casefoam.mkCases(baseCase, caseStructure, caseData,
                 hierarchy='tree', writeDir=baseCaseDir)

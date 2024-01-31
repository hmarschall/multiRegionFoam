import os
import casefoam
import pandas as pd
import matplotlib.pyplot as plt

##################################################
# Retrieve post-processing data
##################################################

# Specify the baseCaseDir and caseStructure as in the genCases.py script
baseCaseDir = 'tests'
caseStructure = [['monolithic', 'partitioned_Aitken'],
                 ['Pr_0.01_Re_500_k_1', 'Pr_0.01_Re_500_k_5', 'Pr_0.01_Re_500_k_20',
                  'Pr_0.01_Re_10000_k_1', 'Pr_0.01_Re_10000_k_5', 'Pr_0.01_Re_10000_k_20',
                  'Pr_100_Re_500_k_1', 'Pr_100_Re_500_k_5', 'Pr_100_Re_500_k_20']]

# Specify a directory in the postProcessing folder
TDir = 'sets/fluid'

# Retrieve data from all tests (specify file name and time folder)
Tdata = casefoam.positional_field(
    TDir, 'DataFile_T.xy', 10, caseStructure, baseCaseDir)

# Rename columns
Tdata.columns = ['x', 'T', 'coupling method', 'parameters']

# Change data type to numeric
numeric_columns = ['x', 'T']
Tdata[numeric_columns] = Tdata[numeric_columns].apply(
    pd.to_numeric, errors='coerce')

# Calculate theta from the temperature
Tdata['theta'] = (Tdata['T'] - 300) / (310 - 300)

##################################################
# Plots
##################################################

# Define a function for common plot settings


def plotCommonSettings():
    plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5),
               borderaxespad=0., prop={'size': 5}, markerscale=1, scatterpoints=1)
    plt.xlabel(r'$x$', fontsize=6),
    plt.ylabel(r'$\theta$', fontsize=6)
    plt.xlim([0, 1])
    plt.xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], [
               '0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'], fontsize=5)
    plt.gcf().set_size_inches(4, 1.5)
    for axis in ['top', 'bottom', 'left', 'right']:
        plt.gca().spines[axis].set_linewidth(0.5)
    plt.gca().tick_params(axis='x', width=0.5)
    plt.gca().tick_params(axis='y', width=0.5)
    plt.savefig(os.getcwd() + "/FOHPResultsPr%sRe%s.png" %
                (Pr, Re), dpi=300, bbox_inches='tight')
    plt.close()


##################################################
# Pr = 0.01 & Re = 500
##################################################

# Data filtering
Pr = 0.01
Re = 500
DATA = Tdata[Tdata['parameters'].str.startswith('Pr_%s_Re_%s_' % (Pr, Re))]
DATApartAitken = DATA[DATA['coupling method'] == 'partitioned_Aitken']
DATAmono = DATA[DATA['coupling method'] == 'monolithic']

# multiRegionFoam results
plt.scatter(DATAmono['x'], DATAmono['theta'], marker='+',
            facecolors="b", linewidth=0.5,  s=10, label='Monolithic')
plt.scatter(DATApartAitken['x'], DATApartAitken['theta'], marker='o', facecolors='None',
            edgecolors="red", linewidth=0.3,  s=7, label='Partitioned-Aitken')

# Analytical and numerical results from Vynnycky et al. (1998)
reference_path = os.getcwd() + '/references/'

ks = ['5', '20']
labels = ['Vynnycky-analytical', '']
for k, label in zip(ks, labels):
    file_name = 'Pr%sRe%sk%sanalytical.csv' % (Pr, Re, k)
    theta = pd.read_csv(reference_path + file_name,
                        delimiter=',', header=None, names=["x", "theta"])
    plt.plot(theta['x'].values, theta['theta'].values,
             'k-', linewidth=0.4, label=label)

ks = ['1', '5', '20']
labels = ['Vynnycky-numerical', '', '']
for k, label in zip(ks, labels):
    file_name = 'Pr%sRe%sk%snumerical.csv' % (Pr, Re, k)
    theta = pd.read_csv(reference_path + file_name,
                        delimiter=',', header=None, names=["x", "theta"])
    plt.plot(theta['x'].values, theta['theta'].values,
             '--', color='grey', linewidth=0.5, label=label)

    # Add k labels
    k_value = int(k)
    plt.text(1.01, theta['theta'].max(),
             f'k={k_value}', color='g', fontsize=5, weight='bold')

plt.ylim([0.2, 1.1])
plt.yticks([0.2, 0.4, 0.6, 0.8, 1], [
           '0.2', '0.4', '0.6', '0.8', '1'], fontsize=5)
plotCommonSettings()


##################################################
# Pr = 100 & Re = 500
##################################################

# Data filtering
Pr = 100
Re = 500
DATA = Tdata[Tdata['parameters'].str.startswith('Pr_%s_Re_%s_' % (Pr, Re))]
DATApartAitken = DATA[DATA['coupling method'] == 'partitioned_Aitken']
DATAmono = DATA[DATA['coupling method'] == 'monolithic']

# multiRegionFoam results
plt.scatter(DATAmono['x'], DATAmono['theta'], marker='+',
            facecolors="b", linewidth=0.5,  s=10, label='Monolithic')
plt.scatter(DATApartAitken['x'], DATApartAitken['theta'], marker='o', facecolors='None',
            edgecolors="red", linewidth=0.3,  s=7, label='Partitioned-Aitken')

# Analytical and numerical results from Vynnycky et al. (1998)
ks = ['1', '5', '20']
labels = ['Vynnycky-analytical', '', '']
for k, label in zip(ks, labels):
    file_name = 'Pr%sRe%sk%sanalytical.csv' % (Pr, Re, k)
    theta = pd.read_csv(reference_path + file_name,
                        delimiter=',', header=None, names=["x", "theta"])
    plt.plot(theta['x'].values, theta['theta'].values,
             'k-', linewidth=0.4, label=label)

ks = ['1', '5', '20']
labels = ['Vynnycky-numerical', '', '']
for k, label in zip(ks, labels):
    file_name = 'Pr%sRe%sk%snumerical.csv' % (Pr, Re, k)
    theta = pd.read_csv(reference_path + file_name,
                        delimiter=',', header=None, names=["x", "theta"])
    plt.plot(theta['x'].values, theta['theta'].values,
             '--', color='grey', linewidth=0.5, label=label)

    # Add k labels
    k_value = int(k)
    plt.text(1.01, theta['theta'].max(),
             f'k={k_value}', color='g', fontsize=5, weight='bold')

plt.ylim([0, 1])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], [
           '0', '0.2', '0.4', '0.6', '0.8', '1'], fontsize=5)
plotCommonSettings()

##################################################
# Pr = 0.01 & Re = 10000
##################################################

# Data filtering
Pr = 0.01
Re = 10000
DATA = Tdata[Tdata['parameters'].str.startswith('Pr_%s_Re_%s_' % (Pr, Re))]
DATApartAitken = DATA[DATA['coupling method'] == 'partitioned_Aitken']
DATAmono = DATA[DATA['coupling method'] == 'monolithic']

# multiRegionFoam results
plt.scatter(DATAmono['x'], DATAmono['theta'], marker='+',
            facecolors="b", linewidth=0.5,  s=10, label='Monolithic')
plt.scatter(DATApartAitken['x'], DATApartAitken['theta'], marker='o', facecolors='None',
            edgecolors="red", linewidth=0.3,  s=7, label='Partitioned-Aitken')

# Analytical and numerical results from Vynnycky et al. (1998)
ks = ['20']
labels = ['Vynnycky-analytical']
for k, label in zip(ks, labels):
    file_name = 'Pr%sRe%sk%sanalytical.csv' % (Pr, Re, k)
    theta = pd.read_csv(reference_path + file_name,
                        delimiter=',', header=None, names=["x", "theta"])
    plt.plot(theta['x'].values, theta['theta'].values,
             'k-', linewidth=0.4, label=label)

ks = ['1', '5', '20']
labels = ['Vynnycky-numerical', '', '']
for k, label in zip(ks, labels):
    file_name = 'Pr%sRe%sk%snumerical.csv' % (Pr, Re, k)
    theta = pd.read_csv(reference_path + file_name,
                        delimiter=',', header=None, names=["x", "theta"])
    plt.plot(theta['x'].values, theta['theta'].values,
             '--', color='grey', linewidth=0.5, label=label)

    # Add k labels
    k_value = int(k)
    plt.text(1.01, theta['theta'].max(),
             f'k={k_value}', color='g', fontsize=5, weight='bold')

plt.ylim([0, 1])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], [
           '0', '0.2', '0.4', '0.6', '0.8', '1'], fontsize=5)
plotCommonSettings()

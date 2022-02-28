"""
Code for generating plots.
"""
import pandas as pd
import os
import matplotlib
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

sns.set_style("darkgrid")


def plot_stairway_Ne_estimate(infile, outfile):
    """
    figure of N(t) for single run of stairwayplot
    """
    nt = pd.read_csv(infile, sep="\t", skiprows=5)
    nt = nt[nt['year'] > 10]
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    ax.plot(nt['year'], nt['Ne_median'], c="red")
    ax.plot(nt['year'], nt['Ne_2.5%'], c='grey')
    ax.plot(nt['year'], nt['Ne_97.5%'], c='grey')
    f.savefig(outfile, bbox_inches='tight')


def plot_compound_Ne_estimate(infiles, outfile):
    """
    figure of N(t) for multiple runs of stairwayplot
    """
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    for infile in infiles:
        nt = pd.read_csv(infile, sep="\t", skiprows=5)
        nt = nt[nt['year'] > 10]
        ax.plot(nt['year'], nt['Ne_median'], c="red")
    f.savefig(outfile, bbox_inches='tight')


def plot_compound_smcpp(infiles, outfile, model, n_samp, generation_time, pop_id=0):
    coal_rate, P, steps = popn_coal_rate(model, pop_id, n_samp, generation_time)
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    for infile in infiles:
        nt = pd.read_csv(infile, usecols=[1, 2], skiprows=0)
        ax.plot(nt['x'], nt['y'], c="red")
    ax.plot(steps, 1/(2*coal_rate), c="black")
    f.savefig(outfile, bbox_inches='tight')


def plot_compound_msmc(infiles, outfile, model):
    f, ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log")
    dfs = []
    for infile in infiles:
        cat = os.path.basename(infile)[0]
        nt = pd.read_csv(infile, usecols=[1, 2], skiprows=0)
        nt['Sample'] = cat
        dfs.append(nt)
    df = pd.concat(dfs, ignore_index=True)
    for label, grp in df.groupby('Sample'):
        grp.plot(x = 'x', y = 'y',ax = ax, label = label)
    f.savefig(outfile, bbox_inches='tight')


def plot_all_ne_estimates(sp_infiles, msmc_infiles, outfile,
                          model, n_samp, generation_time, species,
                          pop_id=0, steps=None): # smcpp_infiles,
    #ddb = model.get_demography_debugger()
    ddb = model.model.debug()
    if steps is None:
        if(len(ddb.epochs) == 1):
            end_time = 100000
        else:
            end_time = int(ddb.epochs[-2].end_time) + 10000
        steps = np.linspace(1, end_time, 1000)
    coal_rate, P = ddb.coalescence_rate_trajectory(steps=steps,
                                                   lineages=n_samp, # changed to lineages
                                                   double_step_validation=False)
    mm = [x for x in ddb.demography.events if "MassMigration" in str(type(x))]
            
    census_size = ddb.population_size_trajectory(steps=steps)
    for m in reversed(mm):
        if m.proportion == 1.0:
            n = (steps > m.time)
            census_size[n, m.source] = census_size[n, m.dest]
        else:
            print("Error: census size estimate requires that MassMigration proportion == 1.0")
            sys.exit(1)
    steps = steps * generation_time
    num_msmc = set([os.path.basename(infile).split(".")[0] for infile in msmc_infiles])
    num_msmc = sorted([int(x) for x in num_msmc])
    f, ax = plt.subplots(1, 2+len(num_msmc), sharex=True, sharey=True, figsize=(14, 7))

    outLines = []
    # for i, infile in enumerate(smcpp_infiles):
    #     nt = pd.read_csv(infile, usecols=[1, 2], skiprows=0)
    #     line1, = ax[0].plot(nt['x'], nt['y'], alpha=0.8)
    #     for j in range(len(nt["x"])):
    #         outLines.append([nt["x"][j], nt["y"][j], "smcpp", "r" + str(i + 1)])
    # ax[0].plot(steps, 1/(2*coal_rate), c="black")
    # ax[0].set_title("smc++")

    for i, infile in enumerate(sp_infiles):
        nt = pd.read_csv(infile, sep="\t", skiprows=5)
        line2, = ax[1].plot(nt['year'], nt['Ne_median'], alpha=0.8)
        for j in range(0, len(nt["year"]), 2):
            outLines.append([nt["year"][j], nt["Ne_median"][j], "sp", "r" + str(i + 1)])

    ax[1].plot(steps, 1/(2*coal_rate), c="black")
    ax[1].set_title("stairwayplot")
    for i, sample_size in enumerate(num_msmc):
        ct = 0
        for infile in msmc_infiles:
            fn = os.path.basename(infile)
            samp = fn.split(".")[0]
            if(int(samp) == sample_size):
                nt = pd.read_csv(infile, usecols=[1, 2], skiprows=0)
                line3, = ax[2+i].plot(nt['x'], nt['y'], alpha=0.8)
                for j in range(len(nt["x"])):
                    outLines.append([nt["x"][j], nt["y"][j], "msmc_" +
                                     str(sample_size), "r" + str(ct + 1)])
                ct += 1
        ax[2+i].plot(steps, 1/(2*coal_rate), c="black")
        ax[2+i].set_title(f"msmc, ({sample_size} samples)")
    plt.suptitle(f"{species}, population id {pop_id}", fontsize=16)
    for i in range(2+len(num_msmc)):
        ax[i].set(xscale="log")
        ax[i].set_xlabel("time (years ago)")
    red_patch = mpatches.Patch(color='black', label='Coalescence rate derived Ne')
    ax[0].legend(frameon=False, fontsize=10, handles=[red_patch])
    ax[0].set_ylabel("population size")
    f.savefig(outfile, bbox_inches='tight', alpha=0.8)

    txtOUT = os.path.join(os.path.dirname(outfile),"_".join([species,model.id,"pop"+str(pop_id),"sizes"])+".txt")
    with open(txtOUT, "w") as fOUT:
        fOUT.write("\t".join([str(x) for x in ["x", "y", "method", "rep"]]) + "\n")
        for i in range(len(steps)):
            fOUT.write("\t".join([str(x) for x in [steps[i],
                                                   1/(2*coal_rate[i]),
                                                   "coal",
                                                   "r1"]]) + "\n")
            fOUT.write("\t".join([str(x) for x in [steps[i],
                                                   census_size[i][pop_id],
                                                   "census",
                                                   "r1"]]) + "\n")
        for i in range(len(outLines)):
            fOUT.write("\t".join([str(x) for x in outLines[i]])+"\n")


def plot_stairwayplot_coalrate(sp_infiles, outfile,
                               model, n_samp, generation_time, species,
                               pop_id=0, steps=None):  # JRA

    #ddb = model.get_demography_debugger()
    ddb = model.model.debug()
    print("This is ddb")
    print(ddb)
    if steps is None:
        end_time = ddb.epochs[-2].end_time + 10000
        steps = np.linspace(1, end_time, end_time+1)
    num_samples = [0 for _ in range(ddb.num_populations)]
    num_samples[pop_id] = n_samp
    print("These are the samples")
    print(num_samples)
    # Change here to extract the populations
    num_samples = {'YRI':20, 'CEU':0, 'CHB':0, 'Neanderthal':0, 'ArchaicAFR':0}
    pop_id = 'YRI'

    coal_rate, P = ddb.coalescence_rate_trajectory(steps=steps,
                                                   lineages=num_samples, # changed from num_samples to lineages 
                                                   double_step_validation=False)
    steps = steps * generation_time
    f, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(7, 7))
    ax.plot(steps, 1/(2*coal_rate), c="black")
    for infile in sp_infiles:
        nt = pd.read_csv(infile, sep="\t", skiprows=5)
        line2, = ax.plot(nt['year'], nt['Ne_median'], alpha=0.8)
    ax.plot(steps, 1/(2*coal_rate), c="black")
    ax.set_title(f"stairwayplot")
    plt.suptitle(f"{species}, population id {pop_id}", fontsize=16)
    ax.set(xscale="log", yscale="log")
    ax.set_xlabel("time (years ago)")
    red_patch = mpatches.Patch(color='black', label='Coalescence rate derived Ne')
    ax.legend(frameon=False, fontsize=10, handles=[red_patch])
    ax.set_ylabel("population size")
    f.savefig(outfile, bbox_inches='tight', alpha=0.8)


def popn_coal_rate(model, pop_id, n_samp, generation_time, steps=None):
    """
    returns tuple (coal_rate, P, steps) for pop_id
    conditional on the model and configuration
    """
    ddb = model.get_demography_debugger()
    if steps is None:
        end_time = ddb.epochs[-2].end_time + 10000
        steps = np.linspace(1, end_time, end_time+1)
    num_samples = [0 for _ in range(ddb.num_populations)]
    num_samples[pop_id] = n_samp
    coal_rate, P = ddb.coalescence_rate_trajectory(steps=steps,
                                                   lineages=num_samples, # changed to lineages
                                                   double_step_validation=False)
    steps = steps * generation_time
    return coal_rate, P, steps


def plot_sfs(s, outfile):
    """
    Plot the SFS for this simulation
    """
    bins = [n + 1 for n in range(len(s[0]))]
    vals = []
    for i in range(len(s)):
        vals.append([int(x) for x in s[i]])
    if len(s) == 2:
        f, ax = plt.subplots(1, 2, sharey=True, tight_layout=True, figsize=(8, 3))
        ax[0].bar(bins, vals[0])
        ax[1].bar(bins, vals[1])
        ax[0].set_title("all sites")
        ax[1].set_title("masked neutral sites")
        ax[0].set_ylabel("counts")
        ax[0].set_xlabel("derived allele frequency")
        ax[1].set_xlabel("derived allele frequency")
        f.savefig(outfile, bbox_inches='tight')
        plt.close()
    else:
        f, ax = plt.subplots(figsize=(3, 3))
        ax.bar(bins, vals[0])
        ax.set_title("all sites")
        ax.set_ylabel("counts")
        ax.set_xlabel("derived allele frequency")
        f.savefig(outfile, bbox_inches='tight')
        plt.close()

def plot_all_dfe_results(mu, seq_len, input, output):
    """
    plot the comparison of dadi and polydfe results for a 
    single model
    """
    #dadi_pop0 = pd.read_csv(input.dadi_res[0], sep="\t")
    #dadi_pop1 = pd.read_csv(input.dadi_res[1], sep="\t")
    #dadi_pop2 = pd.read_csv(input.dadi_res[2], sep="\t")
    polydfe_pop0 = pd.read_csv(input.polydfe_res[0], sep="\t")
    polydfe_pop1 = pd.read_csv(input.polydfe_res[1], sep="\t")
    polydfe_pop2 = pd.read_csv(input.polydfe_res[2], sep="\t")
    #dadi_pop0_Nes = dadi_pop0['theta']/(4*mu*seq_len*0.7)
    #dadi_pop1_Nes = dadi_pop1['theta']/(4*mu*seq_len*0.7)
    #dadi_pop2_Nes = dadi_pop2['theta']/(4*mu*seq_len*0.7)
    #dadi_pop0_Es = dadi_pop0['shape']*dadi_pop0['scale']/(2*dadi_pop0_Nes)
    #dadi_pop1_Es = dadi_pop1['shape']*dadi_pop1['scale']/(2*dadi_pop1_Nes)
    #dadi_pop2_Es = dadi_pop2['shape']*dadi_pop2['scale']/(2*dadi_pop2_Nes)

    ## theta_bar = 4*Ne*mu
    ## S_d = 4*Ne*E(s)
    ## E(s) = mu*S_d/theta_bar
    ## See Fig. 1 in doi: 10.1534/genetics.117.300323
    polydfe_pop0_Es = abs(mu*polydfe_pop0['S_d']/polydfe_pop0['theta_bar'])
    polydfe_pop1_Es = abs(mu*polydfe_pop1['S_d']/polydfe_pop1['theta_bar'])
    polydfe_pop2_Es = abs(mu*polydfe_pop2['S_d']/polydfe_pop2['theta_bar'])

    fig = plt.figure(figsize=(6,6), dpi=300)

    #plt.subplot(2,2,1)
    #plt.ylabel('Shape parameter')
    #plt.xlabel('Population')
    #plt.boxplot([dadi_pop0['shape'], dadi_pop1['shape'], dadi_pop2['shape']])
    #plt.plot(plt.xlim(), [0.19,0.19], c='red')
    #plt.xticks([1,2,3],['YRI','CEU','CHB'])
    #plt.ylim([0.1,0.25])
    #plt.title('dadi')

    #plt.subplot(2,2,2)
    #plt.ylabel('|E(s)|')
    #plt.xlabel('Population')
    #plt.boxplot([dadi_pop0_Es, dadi_pop1_Es, dadi_pop2_Es])
    #plt.plot(plt.xlim(), [0.014,0.014], c='red')
    #plt.xticks([1,2,3],['YRI','CEU','CHB'])
    #plt.ylim([0,0.04])
    #plt.title('dadi')

    plt.subplot(2,2,3)
    plt.ylabel('Shape parameter')
    plt.xlabel('Population')
    plt.boxplot([polydfe_pop0['b'], polydfe_pop1['b'], polydfe_pop2['b']])
    plt.plot(plt.xlim(), [0.19,0.19], c='red')
    plt.xticks([1,2,3],['YRI','CEU','CHB'])
    plt.ylim([0.1,0.25])
    plt.title('polyDFE')

    plt.subplot(2,2,4)
    plt.ylabel('|E(s)|')
    plt.xlabel('Population')
    plt.boxplot([polydfe_pop0_Es, polydfe_pop1_Es, polydfe_pop2_Es])
    plt.plot(plt.xlim(), [0.014,0.014], c='red')
    plt.xticks([1,2,3],['YRI','CEU','CHB'])
    #plt.ylim([0,0.04])
    plt.title('polyDFE')

    fig.tight_layout()
    plt.savefig(output[0])

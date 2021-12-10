import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use("Agg")

def	plot_all_results(mu, seq_len, input, output):
    """
    plot the comparison of dadi and polydfe results for a 
    single model
    """
    dadi_pop0 = pd.read_csv(input.dadi_res[0], sep="\t")
    dadi_pop1 = pd.read_csv(input.dadi_res[1], sep="\t")
    dadi_pop2 = pd.read_csv(input.dadi_res[2], sep="\t")
    polydfe_pop0 = pd.read_csv(input.polydfe_res[0], sep="\t")
    polydfe_pop1 = pd.read_csv(input.polydfe_res[1], sep="\t")
    polydfe_pop2 = pd.read_csv(input.polydfe_res[2], sep="\t")
    dadi_pop0_Nes = dadi_pop0['theta']/(4*mu*seq_len*0.7)
    dadi_pop1_Nes = dadi_pop1['theta']/(4*mu*seq_len*0.7)
    dadi_pop2_Nes = dadi_pop2['theta']/(4*mu*seq_len*0.7)
    dadi_pop0_Es = dadi_pop0['shape']*dadi_pop0['scale']/(2*dadi_pop0_Nes)
    dadi_pop1_Es = dadi_pop1['shape']*dadi_pop1['scale']/(2*dadi_pop1_Nes)
    dadi_pop2_Es = dadi_pop2['shape']*dadi_pop2['scale']/(2*dadi_pop2_Nes)

    ## theta_bar = 4*Ne*mu
    ## S_d = 4*Ne*E(s)
    ## E(s) = mu*S_d/theta_bar
    ## See Fig. 1 in doi: 10.1534/genetics.117.300323
    polydfe_pop0_Es = abs(mu*polydfe_pop0['S_d']/polydfe_pop0['theta_bar'])
    polydfe_pop1_Es = abs(mu*polydfe_pop1['S_d']/polydfe_pop1['theta_bar'])
    polydfe_pop2_Es = abs(mu*polydfe_pop2['S_d']/polydfe_pop2['theta_bar'])

    fig = plt.figure(figsize=(6,6), dpi=300)

    plt.subplot(2,2,1)
    plt.ylabel('Shape parameter')
    plt.xlabel('Population')
    plt.boxplot([dadi_pop0['shape'], dadi_pop1['shape'], dadi_pop2['shape']])
    plt.plot(plt.xlim(), [0.19,0.19], c='red')
    plt.xticks([1,2,3],['YRI','CEU','CHB'])
    plt.ylim([0.1,0.25])
    plt.title('dadi')

    plt.subplot(2,2,2)
    plt.ylabel('|E(s)|')
    plt.xlabel('Population')
    plt.boxplot([dadi_pop0_Es, dadi_pop1_Es, dadi_pop2_Es])
    plt.plot(plt.xlim(), [0.014,0.014], c='red')
    plt.xticks([1,2,3],['YRI','CEU','CHB'])
    #plt.ylim([0,0.04])
    plt.title('dadi')

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

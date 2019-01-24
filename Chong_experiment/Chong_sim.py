import sys
import TCDS.simulation as sim
import TCDS.plotting as plotting
import numpy as np
import matplotlib.pyplot as plt
import os

# number of simulations
nbr_of_sims = int(sys.argv[1])

def chong_experiment(iteration_number):

    # load the files path
    INI_file=sys.argv[2] #
    output_dir_pre=sys.argv[3]+str(iteration_number) #
    output_dir_post= sys.argv[4]+str(iteration_number) #

    # Editing the params.ini file (removing gyrase)
    # open the file
    config = sim.read_config_file(INI_file)
    # modify the gyrase concentration value
    config.set('SIMULATION','GYRASE_CONC', '0.0')
    # and the simulation time
    config.set('SIMULATION','SIM_TIME', '6000')
    # save it
    with open(INI_file, 'w') as configfile:
        config.write(configfile)

    # starting the simulation
    sim.start_transcribing(INI_file, output_dir_pre)

    # Editing the params.ini file (Adding gyrase)
    # open the file
    config = sim.read_config_file(INI_file)
    # modify the gyrase concentration value
    config.set('SIMULATION','GYRASE_CONC', '0.1')
    # and the simulation time 
    config.set('SIMULATION','SIM_TIME', '3000')
    # save it
    with open(INI_file, 'w') as configfile:
        config.write(configfile)

    # resuming the simulation
    sim.resume_transcription(INI_file, output_dir_pre, output_dir_post)


# repeat the simulation 'i' times
for i in range(nbr_of_sims):
    chong_experiment(i)


# create an empty array which will be filled with init_rate_info
init_rate_info_all = np.empty((0,4500))
mean_sigma_info_all = np.empty((0,4500))
sigma_TSS_info_all = np.empty((0,4500))
nbr_rnaps_hooked_info_all = np.empty((0,4500))

# Read the saved files
for i in range(nbr_of_sims):
    ####### Get information about the initiation rate
    # load the npz files (before adding gyrase)
    init_rate_pre = np.load(sys.argv[3]+str(i)+"/all_res/save_tr_info.npz")['tr_info'][0,1,:]
    # load the npz files (after adding gyrase)
    init_rate_post = np.load(sys.argv[4]+str(i)+"/all_res/save_tr_info.npz")['tr_info'][0,1,:]
    # concatenate each pre with it associated post array
    init_rate_info = np.concatenate([init_rate_pre, init_rate_post])
    ## before gathering the arrays retrieved from pre and post in one big array
    ## we will plot the init_rate of each simulation separetly first
    # make te directories in which the generated graphs (of each simulation) wil be saved
    os.makedirs("sim_plots/init_rate_graphs", exist_ok=True)
    ## plot and save
    fig_init_rate = plt.figure()
    plt.plot(range(0, len(init_rate_info)*int(2), int(2)), init_rate_info)
    plt.title("The initiation rate - Sim%d" %i)
    plt.xlabel("Time (s)")
    plt.ylabel("Initiation rate")
    # Save the figure
    fig_init_rate.savefig('sim_plots/init_rate_graphs/init_rate_sim%d.pdf' %i, transparent=True)

    ### Get information about the mean of sigma
    # load the npz files (before adding gyrase)
    mean_sigma_pre = np.load(sys.argv[3]+str(i)+"/all_res/save_sigma_info.npz")["mean_sig_wholeGenome"]
    # load the npz files (after adding gyrase)
    mean_sigma_post = np.load(sys.argv[4]+str(i)+"/all_res/save_sigma_info.npz")["mean_sig_wholeGenome"]
    # concatenate each pre with it associated post array
    mean_sigma_info = np.concatenate([mean_sigma_pre, mean_sigma_post])
    # make te directories in which the generated graphs (of each simulation) wil be saved
    os.makedirs("sim_plots/mean_sigma_graphs", exist_ok=True)
    ## plot and save
    fig_mean_sigma = plt.figure()
    plt.plot(range(0, len(mean_sigma_info)*int(2), int(2)), mean_sigma_info)
    plt.title("The mean of sigma - Sim%d" %i)
    plt.xlabel("Time (s)")
    plt.ylabel("mean of sigma")
    # Save the figure
    fig_mean_sigma.savefig('sim_plots/mean_sigma_graphs/mean_sigma_sim%d.pdf' %i, transparent=True)


    ### sigma in the promoters (TSS)
    sigma_TSS_pre = np.load(sys.argv[3]+str(i)+"/all_res/save_sigma_info.npz")["Barr_sigma_info"].tolist()
    # extract sigma in the TSS
    sigma_TSS_pre = [col[0] for col in sigma_TSS_pre]
    # the same for post
    sigma_TSS_post = np.load(sys.argv[4]+str(i)+"/all_res/save_sigma_info.npz")["Barr_sigma_info"].tolist()
    sigma_TSS_post = [col[0] for col in sigma_TSS_post]
    # concatenate each pre with it associated post array
    sigma_TSS_info = np.concatenate([sigma_TSS_pre, sigma_TSS_post])
    # make te directories in which the generated graphs (of each simulation) wil be saved
    os.makedirs("sim_plots/sigma_tss_graphs", exist_ok=True)
    ## plot and save
    fig_sigma_TSS = plt.figure()
    plt.plot(range(0, len(sigma_TSS_info)*int(2), int(2)), sigma_TSS_info)
    plt.title("The mean of sigma in the promoter - Sim%d" %i)
    plt.xlabel("Time (s)")
    plt.ylabel("sigma in the promoter")
    # Save the figure
    fig_sigma_TSS.savefig('sim_plots/sigma_tss_graphs/sigma_TSS_sim%d.pdf' %i, transparent=True)


    ### numbre of RNAPs which are transcribing (hooked)
    # load the npz files (before adding gyrase)
    nbr_RNAPs_hooked_pre = np.load(sys.argv[3]+str(i)+"/save_nbr_RNAPs_hooked.npz")['nbr_RNAPs_hooked']
    # load the npz files (after adding gyrase)
    nbr_RNAPs_hooked_post = np.load(sys.argv[4]+str(i)+"/save_nbr_RNAPs_hooked.npz")['nbr_RNAPs_hooked']
    # concatenate each pre with it associated post array
    nbr_RNAPs_hooked_info = np.concatenate([nbr_RNAPs_hooked_pre, nbr_RNAPs_hooked_post])
    # make the directories in which the generated graphs (of each simulation) wil be saved
    os.makedirs("sim_plots/nbr_rnaps_hooked_graphs", exist_ok=True)
    ## plot and save
    fig_nbr_RNAPs_hooked = plt.figure()
    plt.plot(range(0, len(nbr_RNAPs_hooked_info)*int(2), int(2)), nbr_RNAPs_hooked_info)
    plt.title("Number of RNAPs which are transcribing - Sim%d" %i)
    plt.xlabel("Time (s)")
    plt.ylabel("Number of hooked RNAPs")
    # Save the figure
    fig_nbr_RNAPs_hooked.savefig('sim_plots/nbr_rnaps_hooked_graphs/nbr_rnaps_hooked_sim%d.pdf' %i, transparent=True)


    ### gather each info in one big array
    # the initiation rate
    init_rate_info_all = np.vstack((init_rate_info_all, init_rate_info))
    # mean of sigma 
    mean_sigma_info_all = np.vstack((mean_sigma_info_all, mean_sigma_info))
    # sigma in the TSS
    sigma_TSS_info_all = np.vstack((sigma_TSS_info_all, sigma_TSS_info))
    # Number of hooked RNAPs
    nbr_rnaps_hooked_info_all = np.vstack((nbr_rnaps_hooked_info_all, nbr_RNAPs_hooked_info))


### Calculate the means
# mean of init rate
moyenne_init_rate = np.mean(init_rate_info_all, axis=0)
# mean of mean sigma
moyenne_mean_sigma = np.mean(mean_sigma_info_all, axis=0)
# mean of sigma in the promoter
moyenne_sigma_TSS = np.mean(sigma_TSS_info_all, axis=0)
# mean of the number of hooked RNAPs
moyenne_nbr_rnaps_hooked = np.mean(nbr_rnaps_hooked_info_all, axis=0)

########################################
## plot the mean of all the simulations
########################################

# Plot the init_rate graph
fig_all_init_rate = plt.figure()
plt.plot(range(0, len(moyenne_init_rate)*int(2), int(2)), moyenne_init_rate)
plt.title("The mean of init_rate")
plt.xlabel("Time (s)")
plt.ylabel("Initiation rate")
# Save the figure
fig_all_init_rate.savefig('sim_plots/%dsim_init_rate.pdf'%nbr_of_sims, transparent=True)

# Plot the mean_sigma graph
fig_all_mean_sigma = plt.figure()
plt.plot(range(0, len(moyenne_mean_sigma)*int(2), int(2)), moyenne_mean_sigma)
plt.title("The mean of mean_sigma")
plt.xlabel("Time (s)")
plt.ylabel("Mean of sigma")
# Save the figure
fig_all_mean_sigma.savefig('sim_plots/%dsim_mean_sigma.pdf'%nbr_of_sims, transparent=True)

# Plot the mean_sigma_TSS graph
fig_all_sigma_TSS = plt.figure()
plt.plot(range(0, len(moyenne_sigma_TSS)*int(2), int(2)), moyenne_sigma_TSS)
plt.title("The mean of sigma in the promoter")
plt.xlabel("Time (s)")
plt.ylabel("Mean of sigma in the TSS")
# Save the figure
fig_all_sigma_TSS.savefig('sim_plots/%dsim_sigma_TSS.pdf'%nbr_of_sims, transparent=True)

# Plot the nbr_rnaps_hooked graph
fig_all_nbr_rnaps_hooked = plt.figure()
plt.plot(range(0, len(moyenne_nbr_rnaps_hooked)*int(2), int(2)), moyenne_nbr_rnaps_hooked)
plt.title("The mean of the number of hooked RNAPs")
plt.xlabel("Time (s)")
plt.ylabel("nbr of RNAPs hooked")
# Save the figure
fig_all_nbr_rnaps_hooked.savefig('sim_plots/%dsim_nbr_rnaps_hooked.pdf'%nbr_of_sims, transparent=True)


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def legend_without_duplicate_labels(ax):
    """
    Remove duplicated labels in legends
    """
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

def set_general_format():
    plt.rcParams.update({'font.size': 9})


class CurrentPlot(object):
    """
    Create a figure that plots protocol, current and state occupancy
    """
    def __init__(self, model):
        super(CurrentPlot, self).__init__()

        plt.rcParams.update({'font.size': 9})

        self.model = model # BindingKinetics model
        
    def add_plot_current(self, signal, voltage_name, pulse_time,
            normal_signal=None, ylim=None):
        # set up figure
        self.fig = plt.figure(figsize=(5, 4))
        gs = self.fig.add_gridspec(3, 1, height_ratios=[1, 2, 2], hspace=0.1)
        self.axs = [self.fig.add_subplot(gs[i, 0]) for i in range(3)]
        
        # plot graph
        self.axs[0].plot(signal.time(), signal[voltage_name]) #, pulse])
        self.axs[1].plot(signal.time(), signal['ikr.IKr'], label='30nM dofetilide') #, pulse])
        self.model.state_occupancy_plot(self.axs[2], signal) #, legend=False) #, pulse=pulse)
        if normal_signal != None:
            self.axs[1].plot(normal_signal.time(), normal_signal['ikr.IKr'], label='control')
        
        # adjust layout
        # self.axs[1].legend()
        self.axs[0].set_ylabel('Voltage (mV)')
        self.axs[1].set_ylabel('Current (A/F)')
        self.axs[2].set_ylabel('State occupancy')

        self.axs[2].set_xlabel('Time (ms)')
        self.axs[2].set_xlim(0, pulse_time)
        if ylim != None:
            self.axs[1].set_ylim(0, ylim)

        for i in range(len(self.axs) - 1):
            self.axs[i].sharex(self.axs[2])
            self.axs[i].tick_params(labelbottom=False)

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.subplots_adjust(hspace=0)


    def add_plot_current_various(self, signals, drug_conc, pulse_time, cmap=None):
        # set up color
        norm = matplotlib.colors.Normalize(0, len(signals))

        # set up figure
        self.fig = plt.figure(figsize=(5, 3))
        gs = self.fig.add_gridspec(2, 1, height_ratios=[1, 2], hspace=0.1)
        self.axs = [self.fig.add_subplot(gs[i, 0]) for i in range(2)]

        # plot graph
        self.axs[0].plot(signals[0].time(), signals[0]['membrane.V'], zorder=-10)
        for i in range(len(signals)):
            self.axs[1].plot(signals[i].time(), signals[i]['ikr.IKr'],
                    label=str(drug_conc[i]) + ' nM', # color=cmap(norm(i)),
                    zorder=-10)

        self.legend = self.axs[1].legend()

        # adjust layout
        self.axs[0].set_ylabel('Voltage (mV)')
        self.axs[1].set_ylabel('Current (A/F)')

        self.axs[1].set_xlabel('Time (ms)')
        self.axs[1].set_xlim(0, pulse_time)

        for i in range(len(self.axs) - 1):
            self.axs[i].sharex(self.axs[1])
            self.axs[i].tick_params(labelbottom=False)

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.subplots_adjust(hspace=0)


    def add_plot_current_various_compare(self, signals1, signals2, drug_conc,
            pulse_time, model_name, cmap=None):
        # set up color
        norm = matplotlib.colors.Normalize(0, len(signals2))

        # set up figure
        self.fig = plt.figure(figsize=(5, 3))
        gs = self.fig.add_gridspec(3, 1, height_ratios=[1, 2, 2], hspace=0.1)
        self.axs = [self.fig.add_subplot(gs[i, 0]) for i in range(3)]

        # plot graph
        self.axs[0].plot(signals1[0].time(), signals1[0]['membrane.V'], 'k', zorder=-10)
        self.axs[0].plot(signals2[0].time(), signals2[0]['membrane.V'], 'k', zorder=-10)
        for i in range(len(signals1)):
            self.axs[1].plot(signals1[i].time(), signals1[i]['ikr.IKr'],
                    label=str(drug_conc[i]) + ' nM', # color=cmap(norm(i)),
                    zorder=-10)
            self.axs[2].plot(signals2[i].time(), signals2[i]['ikr.IKr'],
                    label=str(drug_conc[i]) + ' nM', # color=cmap(norm(i)),
                    zorder=-10)
        self.legend = self.axs[2].legend()

        # adjust layout
        self.axs[0].set_ylabel('Voltage (mV)')
        self.axs[1].set_ylabel('Current (A/F) \n '+ model_name[0])
        self.axs[2].set_ylabel('Current (A/F) \n '+ model_name[1])

        self.axs[2].set_xlabel('Time (ms)')
        self.axs[2].set_xlim(0, pulse_time)

        for i in range(len(self.axs) - 1):
            self.axs[i].sharex(self.axs[2])
            self.axs[i].tick_params(labelbottom=False)

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.subplots_adjust(hspace=0)


    def add_plot_current_pulses(self, signals, pulse_time, start_pulse=0, end_pulse=None, cmap=None):
        
        # wrong end pulse, should take string before the dot
        if end_pulse == None:
            end_pulse = max([int(i[0]) for i in signals.keys_like('membrane.V')]) 
        
        # set up color
        norm = matplotlib.colors.Normalize(0, end_pulse - start_pulse)

        # set up figure
        self.fig = plt.figure(figsize=(5, 3))
        gs = self.fig.add_gridspec(2, 1, height_ratios=[1, 2], hspace=0.1)
        self.axs = [self.fig.add_subplot(gs[i, 0]) for i in range(2)]

        # plot graph
        self.axs[0].plot(signals.time(), signals['membrane.V', 0], zorder=-10)
        for i in range(end_pulse - start_pulse):
            self.axs[1].plot(signals.time(), signals['ikr.IKr', start_pulse + i],
                    color=cmap(norm(i)), zorder=-10)

        self.legend = self.axs[1].legend()

        # adjust layout
        self.axs[0].set_ylabel('Voltage (mV)')
        self.axs[1].set_ylabel('Current (A/F)')

        self.axs[1].set_xlabel('Time (ms)')
        self.axs[1].set_xlim(0, pulse_time)

        for i in range(len(self.axs) - 1):
            self.axs[i].sharex(self.axs[1])
            self.axs[i].tick_params(labelbottom=False)

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.subplots_adjust(hspace=0)


    def add_plot_AP(self, signal, pulse_time, normal_signal=None):
        # set up figure
        self.fig = plt.figure(figsize=(5, 6))
        gs = self.fig.add_gridspec(4, 1, height_ratios=[1, 2, 2, 2], hspace=0.1)
        self.axs = [self.fig.add_subplot(gs[i, 0]) for i in range(4)]

        # plot graph
        self.axs[0].plot(signal.time(), signal['stimulus.i_stim'])
        self.axs[1].plot(signal.time(), signal['membrane.V'], label='30nM dofetilide')
        self.axs[2].plot(signal.time(), signal['ikr.IKr'], label='30nM dofetilide')
        self.model.state_occupancy_plot(self.axs[3], signal) #, legend=False)

        if normal_signal != None:
            self.axs[1].plot(normal_signal.time(), normal_signal['membrane.V'],
                    label='control')
            self.axs[2].plot(normal_signal.time(), normal_signal['ikr.IKr'],
                    label='control')
        self.axs[1].legend()
        self.axs[2].legend()

        # adjust layout
        self.axs[0].set_ylabel('Current stimulus')
        self.axs[1].set_ylabel('Voltage (mV)')
        self.axs[2].set_ylabel('hERG current')
        self.axs[3].set_ylabel('State occupancy')

        self.axs[3].set_xlabel('Time (ms)')
        self.axs[3].set_xlim(0, pulse_time)
        self.axs[2].set_ylim(0, 0.9)

        for i in range(len(self.axs) - 1):
            self.axs[i].sharex(self.axs[3])
            self.axs[i].tick_params(labelbottom=False)

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.subplots_adjust(hspace=0)


    def add_plot_APs(self, *args, pulse_time):
        # set up figure
        self.fig = plt.figure(figsize=(8, 6))
        gs = self.fig.add_gridspec(4, len(args), height_ratios=[1, 2, 2, 2], hspace=0.1, wspace=0.1)
        axs = [[self.fig.add_subplot(gs[i, j]) for i in range(4)] for j in range(len(args))]

        # plot graph
        count = 0
        for i in range(len(args)):
            sgn = args[i]
            axs[i][0].plot(sgn.time(), sgn['stimulus.i_stim'])
            axs[i][1].plot(sgn.time(), sgn['membrane.V'])
            axs[i][2].plot(sgn.time(), sgn['ikr.IKr'])
            self.model.state_occupancy_plot(axs[i][3], sgn, legend=(False if count != 2 else True))
            axs[i][3].set_xlabel('Time (ms)')
            axs[i][3].set_xlim(0, pulse_time)
            count += 1

        # adjust layout
        axs[0][0].set_ylabel('Current stimulus')
        axs[0][1].set_ylabel('Voltage (mV)')
        axs[0][2].set_ylabel('hERG current')
        axs[0][3].set_ylabel('State occupancy')

        for j in range(len(args)):
            for i in range(len(axs[j]) - 1):
                axs[j][i].sharex(axs[j][3])
                axs[j][i].tick_params(labelbottom=False)

        for j in range(len(args) - 1):
            for i in range(len(axs[j + 1])):
                axs[j + 1][i].sharey(axs[0][i])
                axs[j + 1][i].tick_params(labelleft=False)

        plt.tight_layout(pad=0.4)#, w_pad=0.1, h_pad=1.0)
        plt.subplots_adjust(hspace=0)


    def add_plot_AP_various(self, signals, drug_conc, pulse_time, pulse_num):
        # set up figure
        self.fig = plt.figure(figsize=(6, 5.5))
        gs = self.fig.add_gridspec(3, 1, height_ratios=[1, 2, 2], hspace=0.1)
        self.axs = [self.fig.add_subplot(gs[i, 0]) for i in range(3)]

        # colormap
        cmap = matplotlib.cm.get_cmap('viridis')
        norm = matplotlib.colors.Normalize(0, len(signals))
        
        # plot graph
        stim_key = [x for x in signals[0].keys() if x.endswith('stimulus.i_stim')]
        stim_key.sort()
        voltage_key = [x for x in signals[0].keys() if x.endswith('membrane.V')]
        voltage_key.sort()
        current_key = [x for x in signals[0].keys() if x.endswith('ikr.IKr')]
        current_key.sort()
        for pulse in range(pulse_num):
            self.axs[0].plot(signals[0].time() + pulse * max(signals[0].time()),
                signals[0][stim_key[pulse]], color='k')
            for i in range(len(signals)):
                self.axs[1].plot(signals[i].time() + pulse * max(signals[i].time()),
                    signals[i][voltage_key[pulse]],
                        label=str(drug_conc[i]) + 'nM',
                        color=cmap(norm(i)))
                self.axs[2].plot(signals[i].time() + pulse * max(signals[i].time()),
                    signals[i][current_key[pulse]],
                        label=str(drug_conc[i]) + 'nM',
                        color=cmap(norm(i)))

        for i in range(len(self.axs) - 1):
            self.axs[i].sharex(self.axs[2])
            self.axs[i].tick_params(labelbottom=False)

        plt.tight_layout(pad=0.4)#, w_pad=0.5, h_pad=1.0)
        plt.subplots_adjust(hspace=0)


    def add_plot_2AP_various(self, signals, drug_conc, pulse_time, pulse_num, fig_width=6):
    # Used

        # set up figure
        self.fig = plt.figure(figsize=(fig_width, 5.5))
        gs = self.fig.add_gridspec(3, 1, height_ratios=[1, 2, 2], hspace=0.1)
        self.axs = [self.fig.add_subplot(gs[i, 0]) for i in range(3)]

        # colormap
        cmap = matplotlib.cm.get_cmap('viridis')
        norm = matplotlib.colors.Normalize(0, len(signals))
        
        # plot graph
        for pulse in range(pulse_num):
            self.axs[0].plot(signals[0].time() + pulse * max(signals[0].time()),
                signals[0]['stimulus.i_stim', pulse], color='k')
            for i in range(len(signals)):
                self.axs[1].plot(signals[i].time() + pulse * max(signals[i].time()),
                    signals[i]['membrane.V', pulse],
                        label=str(drug_conc[i]) + 'nM',
                        color=cmap(norm(i)))
                self.axs[2].plot(signals[i].time() + pulse * max(signals[i].time()),
                    signals[i]['ikr.IKr', pulse],
                        label=str(drug_conc[i]) + 'nM',
                        color=cmap(norm(i)))

        self.axs[2].legend()
        legend_without_duplicate_labels(self.axs[2])

        # adjust layout
        self.axs[0].set_ylabel('Current stimulus')
        self.axs[1].set_ylabel('Voltage (mV)')
        self.axs[2].set_ylabel('hERG current')

        self.axs[2].set_xlabel('Time (ms)')
        self.axs[2].set_xlim(0, pulse_time)

        for i in range(len(self.axs) - 1):
            self.axs[i].sharex(self.axs[2])
            self.axs[i].tick_params(labelbottom=False)

        plt.tight_layout(pad=0.4)#, w_pad=0.5, h_pad=1.0)
        plt.subplots_adjust(hspace=0)


    def hERG_compare(self, signals_trapping, signals_conductance,
    # Used

            drug_conc, pulse_time, row=1, col=1):
        # set up figure
        self.fig = plt.figure(figsize=(2 * col, 2 * row))
        gs = self.fig.add_gridspec(row, col, hspace=0.3)
        self.axs = [self.fig.add_subplot(gs[i, j]) for i in range(row) for j in range(col)]

        # plot graph
        for i in range(len(drug_conc)):
            self.axs[i].plot(signals_trapping[i].time(),
                    signals_trapping[i]['ikr.IKr'], label='trapping')
            self.axs[i].plot(signals_conductance[i].time(),
                    signals_conductance[i]['ikr.IKr'], label='w/o trapping')
            self.axs[i].set_title('%.1e nM' % drug_conc[i], fontsize=8)

        self.axs[0].legend()

        # adjust layout
        for i in range(row):
            self.axs[i * col].set_ylabel('Current stimulus')
            for j in range(col - 1):
                self.axs[i * col + j + 1].sharey(self.axs[i])
                self.axs[i * col + j + 1].tick_params(labelleft=False)

        for c in range(col):
            self.axs[c + col].set_xlabel('Time (ms)')
            self.axs[c + col].set_xlim(0, pulse_time)
            self.axs[c].sharex(self.axs[c + col])
            self.axs[c].tick_params(labelbottom=False)

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.subplots_adjust(hspace=0)


    def adjust_ticks(self, ax, pulse_time):
        ax.xaxis.set_major_locator(ticker.MultipleLocator(pulse_time / 5))
        xaxis_label = list(ax.get_xticks())
        xaxis_label = ["%d"%(float(i)/1000) for i in xaxis_label]

        ax.set_xticklabels(xaxis_label)

        return ax
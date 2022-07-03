import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

class FigureStructure(object):
    """
    Create structure of a figure with several reference structure
    """
    def __init__(self, figsize=(5, 4), gridspec=(1, 1),
            height_ratios=None, hspace=0.1, wspace=None):
        super(FigureStructure, self).__init__()

        plt.rcParams.update({'font.size': 9})

        self.fig = plt.figure(figsize=figsize)
        self.gridspec = gridspec

        if height_ratios == None:
            height_ratios = [1] + [2] * (self.gridspec[0] - 1)

        gs = self.fig.add_gridspec(*self.gridspec,
                height_ratios=height_ratios, hspace=hspace, wspace=wspace)
        self.axs = [[self.fig.add_subplot(gs[i, j]) for j in range(self.gridspec[1])] for i in range(self.gridspec[0])]

    def legend_without_duplicate_labels(self, ax):
        """
        Remove duplicated labels in legends
        """
        handles, labels = ax.get_legend_handles_labels()
        unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
            
        return unique

    def sharex(self, xlabels, xlim=None):
        """
        Share x-axis for each column and add labels.
        """
        col = self.gridspec[1]
        row = self.gridspec[0]
        
        if len(xlabels) != col:
            raise ValueError('Number of labels must be equal to number of columns.')

        for j in range(col):
            if xlim != None:
                if len(xlim) != col:
                    raise ValueError('Number of limits must be equal to number of columns.')
                self.axs[row - 1][j].set_xlim(xlim[j])

            for i in range(row - 1):
                self.axs[i][j].sharex(self.axs[row - 1][j])
                self.axs[i][j].tick_params(labelbottom=False)
            self.axs[row - 1][j].set_xlabel(xlabels[j])

    def sharey(self, ylabels, ylim=None):
        """
        Share y-axis for each row and add labels
        """

        col = self.gridspec[1]
        row = self.gridspec[0]

        if len(ylabels) != row:
            raise ValueError('Number of labels must be equal to number of rows.')

        for i in range(row):
            if ylim != None:
                if len(ylim) != row:
                    raise ValueError('Number of limits must be equal to number of rows.')
                self.axs[i][0].set_ylim(ylim[i])

            for j in range(col - 1):
                self.axs[i][j + 1].sharey(self.axs[i][0])
                self.axs[i][j + 1].tick_params(labelleft=False)
            self.axs[i][0].set_ylabel(ylabels[i])

    def adjust_ticks(self, ax, pulse_time):
        ax.xaxis.set_major_locator(ticker.MultipleLocator(pulse_time / 5))
        xaxis_label = list(ax.get_xticks())
        xaxis_label = ["%d"%(float(i)/1000) for i in xaxis_label]

        ax.set_xticklabels(xaxis_label)

    def savefig(self, filename):

        plt.subplots_adjust(hspace=0)
        self.fig.savefig(filename, bbox_inches='tight')
        plt.close()

    def reference_structure(self):
        self.fig_structure = {
            'current_state': {
                'figsize': (5, 4),
                'gridspec': (3, 1),
                'row0': 'membrane.V',
                'row0-fn': 'add_single',
                'row1': 'ikr.IKr',
                'row1-fn': 'add_single',
                'row2-fn': 'state_occupancy_plot',
            },}

class FigurePlot(object):
    """
    Create structure of a figure with several reference structure
    """
    def __init__(self):
        super(FigurePlot, self).__init__()

    def add_single(self, ax, log, key, color=None, label=None, alpha=None):
        
        ax.plot(log.time(), log[key], color=color, label=label, alpha=alpha)

    def add_multiple(self, ax, log, key,
            labels=None, color=None):
        """
        Plot overlapping signals
        """
        if color != None:
            norm = matplotlib.colors.Normalize(0, len(log))

        if labels != None and color != None:
            for i in range(len(log)):
                ax.plot(log[i].time(), log[i][key],
                        label=str(labels[i]), color=color(norm(i)), zorder=-10)
        elif color != None:
            for i in range(len(log)):
                ax.plot(log[i].time(), log[i][key], color=color(norm(i)), zorder=-10)
        elif labels != None:
            for i in range(len(log)):
                ax.plot(log[i].time(), log[i][key], label=str(labels[i]), zorder=-10)
        else:
            for i in range(len(log)):
                ax.plot(log[i].time(), log[i][key], zorder=-10)
        
    def add_continuous(self, ax, log, key, start_pulse=0, end_pulse=None,
            label='', color='k'):
        
        if end_pulse == None:
            num_keys = [x for x in log.keys() if x.endswith(key)]
            end_pulse = len(num_keys)

        for pulse in range(end_pulse - start_pulse):
            ax.plot(log.time() + pulse * max(log.time()), log[key, start_pulse + pulse],
                    label=label, color=color, zorder=-10)

    def add_multiple_continuous(self, ax, log, key, start_pulse=0, end_pulse=None,
            labels=None, cmap=None):
        
        if end_pulse == None:
            num_keys = [x for x in log[0].keys() if x.endswith(key)]
            end_pulse = len(num_keys)

        if cmap != None:
            norm = matplotlib.colors.Normalize(0, len(log))

        if labels != None and cmap != None:
            for pulse in range(end_pulse - start_pulse):
                for i in range(len(log)):
                    ax.plot(log[i].time() + pulse * max(log[i].time()),
                            log[i][key, pulse], label=str(labels[i]),
                            color=cmap(norm(i)), zorder=-10)
        elif cmap != None:
            for pulse in range(end_pulse - start_pulse):
                for i in range(len(log)):
                    ax.plot(log[i].time() + pulse * max(log[i].time()),
                            log[i][key, pulse], color=cmap(norm(i)), zorder=-10)
        elif labels != None:
            for pulse in range(end_pulse - start_pulse):
                for i in range(len(log)):
                    ax.plot(log[i].time() + pulse * max(log[i].time()),
                            log[i][key, pulse], label=str(labels[i]), zorder=-10)
        else:
            for pulse in range(end_pulse - start_pulse):
                for i in range(len(log)):
                    ax.plot(log[i].time() + pulse * max(log[i].time()),
                            log[i][key, pulse], zorder=-10)


    def state_occupancy_plot(self, ax, signal_log, model, pulse=None, legend=True, legend_names=None):

        if pulse == None:
             ax.stackplot(signal_log.time(),
                 *[signal_log[s] for s in model.states() \
                         if str(s.parent()) == 'ikr' and s.name() != 'D'],
                 labels = [s.name() for s in model.states() \
                         if str(s.parent()) == 'ikr' and s.name() != 'D'],
                 zorder=-10)
        else:
            ax.stackplot(signal_log.time(),
                *[signal_log[s, pulse] for s in model.states() \
                        if str(s.parent()) == 'ikr' and s.name() != 'D'],
                labels = [s.name() for s in model.states() \
                        if str(s.parent()) == 'ikr' and s.name() != 'D'],
                zorder=-10)
        
        ax.set_xlabel('Time (ms)')
        
        label_list = []
        for t in ax.get_legend_handles_labels():
            label_list.append(t)
       
        if legend_names == None:
            new_list = label_list[1]
            new_list = ['O*' if x == 'Obound' else x for x in new_list]
            new_list = ['C*' if x == 'Cbound' else x for x in new_list]
            new_list = ['IO*' if x == 'IObound' else x for x in new_list]
        else:
            new_list = legend_names

        if legend:
            ax.legend(ncol=2, handles=label_list[0], labels=new_list,
                    loc="lower right", handlelength=1, columnspacing=1,
                    labelspacing=0.3)
        ax.set_rasterization_zorder(0)
        
    def add_overlapping_pulses(self, ax, log, key, start_pulse, end_pulse=None):
        if end_pulse == None:
            num_keys = [x for x in log.keys() if x.endswith(key)]
            end_pulse = len(num_keys)

        for pulse in range(end_pulse - start_pulse):
            ax.plot(log.time(), log[key, start_pulse + pulse], zorder=-10)

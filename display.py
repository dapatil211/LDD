import re
import pprint
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
pp = pprint.PrettyPrinter(indent=2)
file_prefix="outputs/ex"
betas=[".002", ".005", ".01", ".02", ".05", ".1", ".2"]

#betas=[".002", ".005"]

#num_threads=["1", "2", "4", "8", "16", "36"]
num_threads=["1", "2", "4"]

def extract_values():
    values = {
        'serial': {
            'clusters':[],
            'cut_edges':[],
            'time':[]
        },
        'par': {
        }
    }
    for num_thread in num_threads:
        values['par'][num_thread] = {
            'clusters':[],
            'cut_edges':[],
            'time':[]
        }

    for beta in betas:
        filename = file_prefix + "_serial_" + str(beta)
        with open(filename) as f:
            content = f.read().split('\n\n')[1:]
            clusters = 0.
            cut_edges = 0.
            time = 0.
            for line in content:
                nums = re.findall('(?:\d|\.|(?:e\-))+', line)
                clusters += int(nums[3])
                cut_edges += int(nums[4])
                time += float(nums[5])
            clusters /= len(content)
            cut_edges /= len(content)
            time /= len(content)
            values['serial']['clusters'].append(clusters)
            values['serial']['cut_edges'].append(cut_edges)
            values['serial']['time'].append(time)
        for num_thread in num_threads:
            filename = file_prefix + '_par_' + str(num_thread) + '_' + str(beta)
            with open(filename) as f:
                content = f.read().split('\n\n')[1:]
                clusters = 0.
                cut_edges = 0.
                time = 0.
                for line in content:
                    nums = re.findall('(?:\d|\.|(?:e\-))+', line)
                    clusters += int(nums[4])
                    cut_edges += int(nums[5])
                    time += float(nums[6])
                clusters /= len(content)
                cut_edges /= len(content)
                time /= len(content)
                values['par'][num_thread]['clusters'].append(clusters)
                values['par'][num_thread]['cut_edges'].append(cut_edges)
                values['par'][num_thread]['time'].append(time)
    return values

def make_graph(values):
    fig, ax = plt.subplots(sharex=True)
    serial_clusters = values['serial']['clusters']
    ax.plot(betas, serial_clusters, label='Serial')
    for num_thread in num_threads:
        ax.plot(betas, values['par'][num_thread]['clusters'], label=str(num_thread) + " Thread Parallel")
    ax.set(xlabel='Beta', ylabel='Number of clusters')
    ax.set_ylim([0,20])
    ax.set_xscale('log')
    ax.legend(loc='best', fontsize='medium')
    plt.grid(True)
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.get_yaxis().set_major_formatter(ScalarFormatter())
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    #ax.ticklabel_format(style='plain', axis='both')
    fig.savefig('clusters.png')
    plt.show()

    fig, ax = plt.subplots(sharex=True)
    serial_clusters = values['serial']['cut_edges']
    ax.plot(betas, serial_clusters, label='Serial')
    for num_thread in num_threads:
        ax.plot(betas, values['par'][num_thread]['cut_edges'], label=str(num_thread) + " Thread Parallel")
    ax.set(xlabel='Beta', ylabel='Number of Cut Edges')
    ax.set_xscale('log')
    ax.set_ylim([0,35])
    ax.legend(loc='best', fontsize='medium')
    plt.grid(True)
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.get_yaxis().set_major_formatter(ScalarFormatter())
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    #ax.ticklabel_format(style='plain', axis='both')
    fig.savefig('cut_edges.png')
    plt.show()

    fig, ax = plt.subplots(sharex=True)
    serial_clusters = values['serial']['time']
    ax.plot(betas, serial_clusters, label='Serial')
    for num_thread in num_threads:
        ax.plot(betas, values['par'][num_thread]['time'], label=str(num_thread) + " Thread Parallel")
    ax.set(xlabel='Beta', ylabel='Time (s)')
    # ax.set_ylim([0,1])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best', fontsize='medium')
    plt.grid(True)
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.get_yaxis().set_major_formatter(ScalarFormatter())
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    #ax.ticklabel_format(style='plain', axis='both')
    fig.savefig('time.png')
    plt.show()

def main():
    values = extract_values()
    make_graph(values)

if __name__ == "__main__":
    main()
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def rand_jitter(arr):
    """
    Add random jitter to an array of values.
    Args:
        arr (numpy.ndarray): The input array to which random jitter will be added.
    
    Returns:
        numpy.ndarray: The array with added random jitter.
    """
    stdev = .01 * (max(arr) - np.percentile(arr, 10))
    return arr + np.random.randn(len(arr)) * stdev


def plot_readbqs(read_bqs, title="Base call qualities", filename="readbasequal-waterfall.png"):
    """
    Plot the base call qualities of reads as waterfall plot.
    
    Args:
        read_bqs (list of list): A list of lists containing base call qualities for each read.
        title (str): The title for the plot.
    
    Returns:
        None: The function saves the plot as a PNG file.
    """

    # Define two colors for the sequential colormap
    start_color = "#c63a1e"  # red for low quality
    end_color   = "#50b7bc"  # teal for high quality
    two_color_cmap = LinearSegmentedColormap.from_list("two_color_sequential", [start_color, end_color])
    two_color_cmap.set_bad('#efefef') # a sliglhtly darker grey for nan values (background)

    # Pad all lists in read_bqs to the highest length by filling with nan values
    max_length = max(len(bq) for bq in read_bqs)
    read_bqs = [bq + [float('nan')] * (max_length - len(bq)) for bq in read_bqs]

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 7)
    ax.imshow(read_bqs, aspect='auto', interpolation='nearest', cmap=two_color_cmap)
    ax.set_xlabel('Repeat position', fontsize=14)
    ax.set_ylabel('Read number', fontsize=14)
    ax.set_title(title, fontsize=14)
    plt.colorbar(ax.imshow(read_bqs, aspect='auto', cmap=two_color_cmap))
    plt.savefig(filename, dpi=300)


def plot_readbms(read_bms, title="Base methylation", filename="readbasemeth-waterfall.png"):
    """
    Plot the methylation levels of each base in all the reads as a waterfall plot.
    
    Args:
        read_bqs (list of list): A list of lists containing methylation levels for each read.
        bam_id (str): The ID of the BAM file.
        filename (str): The filename for the saved plot.
    
    Returns:
        None: The function saves the plot as a PNG file.    
    """

    # Define two colors for the sequential colormap
    vcolors = [
        (0.0, "#cecece"),    # grey for bases where methylation is not called
        (0.003, "#50b7bc"),  # teal for no methylation no methylation score starts at 1/256
        (0.5, "white"),        # white for 50% methylation
        (1, "#c63a1e"),      # red  for 100% methylation
    ]
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", vcolors)
    custom_cmap.set_bad('#efefef') # a sliglhtly darker grey for nan values (background)

    # Pad all lists in read_bqs to the highest length by filling with nan values
    max_length = max(len(bm) for bm in read_bms)
    read_bms = [bm + [float('nan')] * (max_length - len(bm)) for bm in read_bms]

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 7)
    ax.imshow(read_bms, aspect='auto', interpolation='nearest', cmap=custom_cmap)
    ax.set_xlabel('Repeat position', fontsize=14)
    ax.set_ylabel('Read number', fontsize=14)
    ax.set_title(f'Methylation', fontsize=14)
    fig.suptitle(title, fontsize=20)
    plt.colorbar(ax.imshow(read_bms, aspect='auto', interpolation='nearest', cmap=custom_cmap))
    plt.savefig(filename, dpi=300)


def plot_bqsvsbms(read_bqs, read_bms, title, filename="readbasequal-methcall-waterfall.png"):
    """
    Plot the base call qualities and methylation levels of each base in all the reads as a waterfall plot.
    
    Args:
        read_bqs (list of list): A list of lists containing base call qualities for each read.
        read_bms (list of list): A list of lists containing methylation levels for each read.
        title (str): The title for the plot.
        filename (str): The filename for the saved plot.
    
    Returns:
        None: The function saves the plot as a PNG file.
    """

    plt.rcParams['font.family'] = 'Arial'

    max_length = max(len(bq) for bq in read_bqs)
    read_bqs = [bq + [float('nan')] * (max_length - len(bq)) for bq in read_bqs]
    read_bms = [bm + [float('nan')] * (max_length - len(bm)) for bm in read_bms]

    fig, ax = plt.subplots(ncols=2, nrows=1)
    fig.set_size_inches(20, 7)

    start_color = "#da5a2a"  # White
    end_color =   "#ffdfb9"    # Blue
    two_color_cmap = LinearSegmentedColormap.from_list("two_color_sequential", [start_color, end_color])
    two_color_cmap.set_bad('#ffffff') #

    ax[0].imshow(read_bqs, aspect='auto', interpolation='nearest', cmap=two_color_cmap)
    ax[0].set_xlabel('Repeat position', fontsize=18, fontname='Arial')
    ax[0].set_ylabel('Read number', fontsize=18, fontname='Arial')
    ax[0].set_title(f'Base call qualities', fontsize=20, fontname='Arial')
    ax[0].tick_params(axis='both', labelsize=14)
    cbar_bq = plt.colorbar(ax[0].imshow(read_bqs, aspect='auto', interpolation='nearest', cmap=two_color_cmap, vmin=0, vmax=40, label='Base quality'))
    cbar_bq.ax.tick_params(labelsize=14)
    cbar_bq.set_label('Base quality', fontsize=14, fontname='Arial')

    vcolors = [
        (0.0, "#eeeeee"),  # White at 0
        (0.001, "#50b7bc"),   # Blue at 1
        (0.5, "white"),   # Blue at 1
        (1, "#c63a1e"),  # Red at 0.001
    ]
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", vcolors)
    custom_cmap.set_bad('#ffffff') 

    # for i, a in enumerate(read_bms): print(i, len(a))
    ax[1].imshow(read_bms, aspect='auto', interpolation='nearest', cmap=custom_cmap)
    ax[1].set_xlabel('Repeat position', fontsize=18, fontname='Arial')
    ax[1].set_ylabel('Read number', fontsize=18, fontname='Arial')
    ax[1].set_title(f'Methylation', fontsize=20, fontname='Arial')
    ax[1].tick_params(axis='both', labelsize=14)
    cbar_bm = plt.colorbar(ax[1].imshow(read_bms, aspect='auto', interpolation='nearest', cmap=custom_cmap, vmin=0, vmax=1, label='Methylation'))
    cbar_bm.ax.tick_params(labelsize=14)
    cbar_bm.set_label('Methylation', fontsize=14, fontname='Arial')
    
    fig.suptitle(f'{title}', fontsize=24, fontname='Arial')
    plt.savefig(filename, dpi=300)
    plt.close()

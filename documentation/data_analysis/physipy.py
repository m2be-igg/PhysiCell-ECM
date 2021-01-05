import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats
from pathlib import Path
from scipy import io as sio
import xml.etree.ElementTree as ET

def get_cell_data(timestep, folder_name, variables='all'):
    """Returns a dictionary with the cell output data for the selected variables.

    Parameters
    ----------
    timestep : int
        The time point at which the output was recorded
    folder_path: Path
        The path to the folder where the output (.mat) files are stored
    variables : list
        The variables to be extracted from the output files. If variables
        are not defined, all the available outputs will be saved.
    """

    # All possible output variables written by PhysiCell
    data_labels = [
        'ID',
        'position_x', 'position_y', 'position_z',
        'total_volume',
        'cell_type',
        'cycle_model', 'current_phase', 'elapsed_time_in_phase',
        'nuclear_volume', 'cytoplasmic_volume',
        'fluid_fraction', 'calcified_fraction',
        'orientation_x', 'orientation_y', 'orientation_z',
        'polarity',
        'migration_speed',
        'motility_vector_x', 'motility_vector_y', 'motility_vector_z',
        'migration_bias',
        'motility_bias_direction_x', 'motility_bias_direction_y', 'motility_bias_direction_z',
        'persistence_time',
        'motility_reserved'
    ]

    # Create path name
    time_str = str(timestep).zfill(8)
    file_name = 'output{}_cells_physicell.mat'.format(time_str)
    path_name = folder_name / file_name

    # Read output file
    cell_data = sio.loadmat(path_name)['cells']

    # Select and save the variables of interest
    variables_indexes = [data_labels.index(var) for var in variables]
    cells = {var: cell_data[index, :]
             for var, index in zip (variables, variables_indexes)}

    return cells


def get_me_data(timestep, substance, base_file):
    """Returns an array with the substance concentrations at the middle plane of the domain.

    Parameters
    ----------
    timestep : int
        The time point at which the output was recorded.
    substance : string
        The substance to be quantified.
    base_file: Path
        The path to the folder where the output (.mat/.xml) files are stored.
    """

    # Create path
    xml_file = base_file / 'output0000000{}.xml'.format(timestep)
    me_file = base_file / 'output0000000{}_microenvironment0.mat'.format(timestep)

    # Open XML file to get list of variables
    tree = ET.parse(xml_file)
    root = tree.getroot()
    me_node = root.find('microenvironment')
    me_node = me_node.find('domain')
    mesh_node = me_node.find('mesh')
    variables_node = me_node.find('variables')
    var_children = variables_node.findall('variable')
    variables = [var.get('name') for var in var_children]

    # Get x, y and z coordinates
    # X coordinates
    coord_str = mesh_node.find('x_coordinates').text
    delimiter = mesh_node.find('x_coordinates').get('delimiter')
    x_coords = np.array(coord_str.split(delimiter), dtype=np.float)
    # Y coordinates
    coord_str = mesh_node.find('y_coordinates').text
    delimiter = mesh_node.find('y_coordinates').get('delimiter')
    y_coords = np.array(coord_str.split(delimiter), dtype=np.float)
    # Z coordinates
    coord_str = mesh_node.find('z_coordinates').text
    delimiter = mesh_node.find('z_coordinates').get('delimiter')
    z_coords = np.array(coord_str.split(delimiter), dtype=np.float)
    z_middle_point = int(len(z_coords) / 2)

    # Define the shape of the output array
    data_shape = (len(x_coords), len(y_coords))

    # Load substance data
    me_data = sio.loadmat(me_file)['multiscale_microenvironment']

    # Select the data corresponding to the chosen substance
    substance_index = variables.index(substance)
    substance_data = me_data[substance_index + 4, me_data[2, :] == z_coords[z_middle_point]]

    # Reshape output array to match the x and y coordinates
    substance_data = np.reshape(substance_data, data_shape)

    return substance_data


def cells_in_area_of_interest(cells, coordinate, size):
    """Selects the cells that are placed in a reduced area.

    Based on the passed coordinate, selects the cells that are placed
    between -size and size according to that axis.

    Parameters
    ---------
    cells : DataFrame
        DataFrame containing the cells.
    coordinate : string
        The coordinate in which to apply the interval.
    size : float
        The size of the interval to be applied (in one direction).
    """

    position = 'position_{}'.format(coordinate)
    cells = cells[cells[position] > -size]
    cells = cells[cells[position] < size]

    return cells


def set_cell_view_labels(ax, index, concentrations):
    """Sets axes labels to predefined settings.


    Parameters
    ---------
    ax : ax
        The ax object to be stylized.
    index: int
        The index of the ax object.
    concentrations: array
        The concentration values to print as titles.
    """

    # Title
    ax.set_title("Collagen concentration: {} [mg/mL]".format(concentrations[index]),
                 y=1.05, fontsize=15)
    # Axes labels
    ax.set_xlabel("X Position [$\mu$m]", labelpad=15, fontsize=15)
    if index == 0:
        ax.set_ylabel("Y Position [$\mu$m]", labelpad=10, fontsize=15)


def set_cell_view_style_axes(ax, axes_lim=400):
    """Sets axes style to predefined settings.

    Includes a light gray dashed grid and all spines set as black,
    with ticks pointing inwards.

    Parameters
    ---------
    ax : ax
        The ax object to be stylized.
    axes_lim : int
        The limits of the x and y axes (assuming it goes from -axes_lim to axes_lim)
    """

    ax.set_ylim(-axes_lim, axes_lim)
    ax.set_xlim(-axes_lim, axes_lim)

    # Grid
    ax.yaxis.grid(color='lightgray', linestyle='--', linewidth=0.5, zorder=-1)
    ax.xaxis.grid(color='lightgray', linestyle='--', linewidth=0.5, zorder=-1)

    # Ticks
    ax.tick_params(axis="y", direction="in", right=True)
    ax.tick_params(axis="x", direction="in", top=True)

    # Spines
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.5)
        spine.set_edgecolor('black')


def draw_cells_as_circles(cells, ax, cell_radius=8):
    """Draws cells with the assumed circular geometry.

    Parameters
    ---------
    cells : DataFrame
        DataFrame containing the cells to be plotted.
    ax : ax
        The ax object in which the cells will be drwan.
    cell_radius: int
        Radius of the cells. If not specified, it is defined as 8 microns.
    """

    # Draw cells with a circular geometry
    for x_pos, y_pos in zip(cells['position_x'], cells['position_y']):
        # Create sphere object
        cell_contour = plt.Circle((x_pos, y_pos),
                                  cell_radius,
                                  facecolor='C0', edgecolor='black')
        # Plot object
        ax.add_patch(cell_contour)

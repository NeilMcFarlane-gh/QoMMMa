#!/usr/bin/python3

"""
// This is a helper QoMMMa Python file, which contains the functions used for a growing string method based calculation. //
"""
# Global imports.
from time import asctime
import numpy as np
import os
import shutil


def gsmlog(s, usrdir):
    """
    
    // Function which writes the progress or error message of the given QoMMMa calculation to the GSM log file QoMMMa8_GSM.log. //
    
    Arguments
    ----------
    s : string
        The message to be written to the file QoMMMa8_GSM.log.
    usrdir : string
        The user directory.

    """
    
    # Simply opens QoMMMa8_GSM.log and adds the string, s.
    fd = open(usrdir + '/QoMMMa8_GSM.log', 'a') 
    fd.write(s)
    fd.write('\n')
    fd.write('\n')
    fd.close()
    
def gsmend(s, usrdir):
    """
    
    // Function which ends the given GSM QoMMMa calculation and writes the error to the log file QoMMMa8_GSM.log. //
    
    Arguments
    ----------
    s : string
        The message to be written to the file QoMMMa8.log.
    usrdir : string
        The user directory.

    """
 
    # The time the program ends is recorded along with an error message.
    hlp = asctime()
    qomlog(s + '\n' + '  QoMMMa job ends at ' + str(hlp), usrdir)
    
    # Ending QoMMMa GSM job.
    sys.exit()
    
def read_driving(usrdir):
    """
    
    // Function which generates the array of driving coordinates from the user-specified file driving_coords.in. //
    
    Arguments
    ----------
    usrdir : string
        The user directory.

    """
    
    # Opening the input file and reading in all lines.
    with open('driving_coords.in', 'r') as file:
        lines = file.read().splitlines()

    # Now, a list of lists is created containing each of the driving coordinates.
    # The first four integers define the coordinate, and the final integer defines the direction (1 for increase, -1 for decrease).
    driving_coords = []
    for i in lines:
        driving_temp = i.split()
        map_object = map(int, driving_temp)
        driving_lst = list(map_object)
        driving_coords.append(driving_lst)

    return driving_coords
    
def DE_get_tangent(frontier_dirs, current_nodes, total_nodes, usrdir):
    """
    
    // Function which generates the primitive internal coordinate tangent which describes the reaction pathway. //
    // This function is only used during the growth-phase of GSM. //
    
    Arguments
    ----------
    frontier_dirs : list
        List of the two 'closest' node directories (called frontier nodes) which are always used to define the tangent in the growth phase.
    current_nodes : integer
        Integer which represents the current nodes along the string - used to define the magnitude of the step.
    total_nodes : integer
        Integer which represents the total number of nodes (including reactant and product nodes) which the string should have.
    usrdir : string
        The user directory.

    """
    
    # First, read in the primitive internal coordinates from the frontier node directories.
    prims_r = []
    prims_p = []
    for dir in frontier_dirs:
        os.chdir(dir)
        with open("prims", 'r') as file:
            n_prims = int(file.readline())
            for i in range(1,n_prims):
                if "nodeR" in dir:
                    prims_r[i] = float(file.readline())
                elif "nodeP" in dir:
                    prims_p[i] = float(file.readline())
                    
    # Now, the tangent can be calculated simply as the difference between these primitive internal coordinates.
    tangent = prims_p - prims_r

    # The magnitude of the stepsize is scaled based on the current number of nodes.
    if (total_nodes - current_nodes) > 1:
        stepsize = 1./float(total_nodes - current_nodes)
    else:
        stepsize = 0.5
    tangent *= stepsize
        
    # If a given primitive is below a certain threshold, then it is set to zero and hence not constrained in the Fortran optimisation.
    for i in tangent:
        if (i < 1E-4):
            tangent[i] = 0.0  
    # to-do: modify tangent so that it is of the correct length.
    # it should only be values which are not zero, and the primitive coordinate definitions corresponding to these terms should be retained.
        
    return tangent
    
def DE_add_nodes(frontier_dirs, new_frontier_dirs, tangent, usrdir, to_add=2):
    """
    
    // Function which creates two (or one) new qommma.in files which essentially represent two (or one) new frontier nodes. //
    // The geometries from the previous two frontier nodes are also copied over to the new directories. //
    
    Arguments
    ----------
    frontier_dirs : list
        List of the two current frontier nodes which are about to be replaced.
    new_frontier_dirs : list
        List of the two (or one) new frontier node directories which are about to be added.
    tangent : list
        The primitive internal coordinate tangent which is used to define the geometry and constraints applied to the new nodes.
    usrdir : string
        The user directory.
    to_add : kwarg
        Integer which represents whether 2 (normal node addition) or 1 (final central node) are to be added to the string.

    """
    
    # First, copy over the qommma.in and geometry files from the previous frontier nodes to the new frontier nodes.
    source = frontier_dirs[1]
    inpf_s = source + '/qommma.in'
    geom_s = source + '/geom.xyz'
    destination = new_frontier_dirs[1]
    inpf_d = destination + '/qommma.in'
    geom_d = destination + '/geom.xyz'  
    shutil.copy(source, destination)        
    if (to_add == 2):
        source = frontier_dirs[2]
        inpf_s = source + '/qommma.in'
        geom_s = source + '/geom.xyz'
        destination = new_frontier_dirs[2]
        inpf_d = destination + '/qommma.in'
        geom_d = destination + '/geom.xyz'  
        shutil.copy(source, destination)
    
    # Now, the new qommma.in files in the new directories are modified to include the new tangent.
    # to-do.....
            
    
def DE_reparam_g(node_dirs, current_nodes, total_nodes, usrdir):   
    """
    
    // Function which reparameterises the string during the growth phase. //
    // This ensures even spacing of the nodes along the tangent. //
    
    Arguments
    ----------
    node_dirs : list
        List of all the node directories which have been added to the string so far.
    current_nodes : integer
        Integer which represents the current nodes along the string..
    total_nodes : integer
        Integer which represents the total number of nodes (including reactant and product nodes) which the string should have.
    usrdir : string
        The user directory.

    """
    
    #returns nothing, but makes new files for the reparameterisation and calls fortran code to update the geometry.
    
def SE_get_tangent(frontier_dir, driving_coords, usrdir):
    """
    
    // Function which generates the primitive internal coordinate driving coordinate which describes the reaction pathway. //
    // This function is only used during the growth-phase of SE-GSM. //
    // The magnitude of the tangent is scaled based essentially on trial-and-error, but works for most cases. //
    
    Arguments
    ----------
    frontier_dir : string
        String containing the node (called frontier node) 'closest' to the stationary point.
    driving_coords : list
        The primitive internal coordinates which the user has specified to be likely involved in the reaction pathway.
    usrdir : string
        The user directory.

    """
    
    # Iterate through each driving coordinate and scale each driving coordinate.
    # The step sizes for these driving coordinates are essentially arbitrary.
    # They are made large enough that the step is reasonable, but not too large that the difference between nodes is too great.
    # In driving_coords, the first four integers define the coordinate, and the final integer defines the direction (1 for increase, -1 for decrease).
    tangent = []
    for indice in range(1,len(driving_coords)):
        i = driving_coords[indice]
        driving = i[0:3]
        direction = i[4]
        
        # Now, find out what type of primitive internal coordinate each driving coordinate is by counting the number of zeroes.
        zero_count = 0
        for j in driving:
            if (j == 0):
                zero_count += 1
        
        if (zero_count == 2): # Bond
            tangent[indice] = 0.1 * direction # Angstroms...
        elif (zero_count == 1): # Angle
            tangent[indice] = 0.0872665 * direction # Radians...
        elif (zero_count == 0): # Dihedral torsion
            tangent[indice] = 0.0872665 * direction # Radians...
           
    return tangent

def SE_get_final_tangent(nodeR_dir, frontier_dir, driving_coords, usrdir):
    """
    
    // Function which generates the primitive internal coordinate driving coordinate which describes the reaction pathway. //
    // This function is only used during the growth-phase of SE-GSM, and specifically for the final 2 - 4 nodes. //
    // The magnitude of the tangent is scaled based on how much of a difference is found between the frontier node and the reactant node. //
    // This is approximate, and assumes that the energy profile is roughly symmetrical (not usually true), but works for growing the string. //
    
    Arguments
    ----------
    frontier_dir : string
        String containing the initial reactant node directory.
    frontier_dir : string
        String containing the node (called frontier node) 'closest' to the stationary point.
    driving_coords : list
        The primitive internal coordinates which the user has specified to be likely involved in the reaction pathway.
    usrdir : string
        The user directory.

    """
    
    return tangent
    
def SE_add_node(frontier_dir, new_frontier_dir, tangent, driving_coords, usrdir):
    """
    
    // Function which creates a new qommma.in files which essentially represents the new frontier node. //
    // The geometry from the previous frontier node is also copied over to the new directory. //
    
    Arguments
    ----------
    frontier_dir : string
        String of the current frontier node which is about to be replaced.
    new_frontier_dir : string
        String of the new frontier node directory which is about to be added.
    tangent : list
        The primitive internal coordinate tangent which is used to define the geometry and constraints applied to the new node.
    driving_coords : list
        The primitive internal coordinates which the user has specified to be likely involved in the reaction pathway.
    usrdir : string
        The user directory.

    """
    
    # First, copy over the qommma.in and geometry file from the previous frontier node to the new frontier node.
    source = frontier_dir
    inpf_s = source + '/qommma.in'
    geom_s = source + '/geom.xyz'
    destination = new_frontier_dir
    inpf_d = destination + '/qommma.in'
    geom_d = destination + '/geom.xyz'  
    shutil.copy(source, destination)  

    # Now, prepare the new input prim_constrain_lst for qommma.in using the tangent and the primitive driving coordinates.
    prim_constrain_lst = []
    for indice in range(1,len(driving_coords)):
        i = driving_coords[indice]
        driving = i[0:3]
        dq = tangent[indice]
        to_append = [dq] + driving
        prim_constrain_lst.append(to_append)
    
    # Now, the new qommma.in file in the new directory is modified to include the new tangent.
    with open(inpf_d, 'r') as file:
        inp_data = file.readlines()
    line_num = 0
    for line in inp_data:
        line_num += 1
        if "prim_constrain_lst" in line:
            inp_data[line_num] = prim_constrain_lst
    with open(inpf_d, 'w') as file:
        file.writelines(inp_data)
        

def SE_add_final_nodes(frontier_dir, new_frontier_dirs, tangent, usrdir):
    """
    
    // Function which creates a news qommma.in files for the final nodes which essentially represent the new frontier nodes. //
    // The geometry from the previous frontier node is used as the starting point for all the new nodes and also copied over to the new directories. //
    
    Arguments
    ----------
    frontier_dir : string
        String of the current frontier node which is about to be replaced.
    new_frontier_dirs : list
        List of the 2 - 4 new frontier node directories which are about to be added.
    tangent : list
        The primitive internal coordinate tangent which is used to define the geometry and constraints applied to the new nodes.
    usrdir : string
        The user directory.

    """
    
    #returns nothing, but makes and moves files.
    
def SE_check_delE():
    # not sure, yet...
    return
    
def get_tangents_opt(node_dirs, usrdir):
    """
    
    // Function which generates the set of primitive internal coordinate tangents which describe the reaction pathway. //
    // This function is only used during the optimisation-phase of GSM, and for both the single- and double-ended variant. //
    
    Arguments
    ----------
    node_dirs : list
        List of all node directories which are used to define the tangents.
    usrdir : string
        The user directory.

    """
    
    return tangent_list
    
def gen_input_opt(node_dirs, tangent_list, usrdir):
    """
    
    // Function which creates a series of new qommma.in files for all the nodes. //
    // This function is only used during the optimisation-phase of GSM, and for both the single- and double-ended variant. //
    
    Arguments
    ----------
    node_dirs : list
        List of all node directories which are used to define the tangents.
    tangent_list : array
        Array of all the primitive internal coordinate tangents which connect the nodes and create the reaction pathway.
    usrdir : string
        The user directory.

    """
    
    #returns nothing, but makes and moves files.
    

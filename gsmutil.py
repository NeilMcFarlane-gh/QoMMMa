#!/usr/bin/python3

"""
// This is a helper QoMMMa Python file, which contains the functions used for a growing string method based calculation. //
"""
# Global imports.
from time import asctime


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
    
    #returns nothing, but makes and moves files.
    
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
    
def SE_add_node(frontier_dir, new_frontier_dir, tangent, usrdir):
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
    usrdir : string
        The user directory.

    """
    
    #returns nothing, but makes and moves files.

def SE_add_final_nodes(frontier_dir, new_frontier_dirs, tangent, usrdir):
    """
    
    // Function which creates a new qommma.in files which essentially represents the new frontier node. //
    // The geometry from the previous frontier node is also copied over to the new directory. //
    
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
    

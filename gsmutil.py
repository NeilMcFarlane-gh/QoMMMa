#!/usr/bin/python3

"""
// This is a helper QoMMMa Python file, which contains the functions used for a growing string method based calculation. //
"""
# Global imports.
from time import asctime
import numpy as np
import os
import shutil
import sys
import math

# Local imports.
import qommma

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
    gsmlog(s + '\n' + '  QoMMMa GSM job ends at ' + str(hlp), usrdir)
    
    # Ending QoMMMa GSM job.
    sys.exit()

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
    # These should have already been checked to be the same definition so maths can be performed on them.
    prims_i = []
    prim_list_i = []
    prims_j = []
    prim_list_j = []
    for indice,dir in enumerate(frontier_dirs):
        os.chdir(dir)
        with open("prims", 'r') as pr, open("prim_list", 'r') as pr_ls:
            n_prims = int(pr.readline())
            pr_ls.readline()
            for i in range(1,n_prims):
                if (indice == 1):
                    prims_i.append(float(pr.readline()))
                    prim_list_i.append(pr_ls.readline())
                elif (indice == 2):
                    prims_j.append(float(pr.readline()))
                    prim_list_j.append(pr_ls.readline())
                    
    # Now, the tangent can be calculated simply as the difference between these primitive internal coordinates.
    tangent_temp = []
    for p1, p2 in zip(prims_i, prims_j):
        p = p1 - p2
        tangent_temp.append(p)

    # The magnitude of the stepsize is scaled based on the current number of nodes.
    if (total_nodes - current_nodes) > 1:
        stepsize = 1./float(total_nodes - current_nodes)
    else:
        stepsize = 0.5 # for the central node, if using an odd number of nodes.
    tangent_temp *= stepsize
        
    # If a given primitive is below a certain threshold, then it is set to zero and hence not constrained in the Fortran optimisation.
    # If a given primitive is below a certain threshold, then it is set to zero and hence not constrained in the Fortran optimisation.
    for i in range(0,n_prims):
        # Store the tangent and primitive associated with it.
        temp_tan = tangent_temp[i]
        temp_prim = prim_list_i[i].split()
                
        # Store the type of primitive internal coordinate.
        zero_count = 0
        for prim in temp_prim:
            if int(prim) == 0:
                zero_count += 1

        # Setting appropriate element to zero.        
        if (abs(i) < 1E-2):
            tangent_temp[i] = 0.0  

        ################################################################################
        # TEMPORARY TEST FIX - means that valence and dihedral angles are not included.#
        ################################################################################
        if (zero_count == 0) or (zero_count == 1):
            tangent_temp[i] = 0.0
            
     # We should remove components of the tangent which are zero and save the relevant primitive definitions corresponding to the tangents.
    tangent = []
    tangent_prims = []
    for counter,tang in enumerate(tangent_temp):
        if tang > 0.0:
            tangent.append(tang)
            tangent_prims.append(prim_list_r(counter))
        
    return tangent, tangent_prims
    
def DE_add_nodes(frontier_dirs, new_frontier_dirs, tangent, tangent_prims, usrdir, to_add=2):
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
    tangent_prims : list
        The primitive internal coordinate definititions for the tangent which is used to define the geometry and constraints applied to the new nodes.
    usrdir : string
        The user directory.
    to_add : kwarg
        Integer which represents whether 2 (normal node addition) or 1 (final central node) are to be added to the string.

    """

    # First, copy over the qommma.in and geometry files from the previous frontier nodes to the new frontier nodes.
    source = frontier_dirs[1]
    inpf_s_i = source + '/qommma.in'
    geom_s_i = source + '/geom1.xyz'
    destination = new_frontier_dirs[1]
    inpf_d_i = destination + '/qommma.in'
    geom_d_i = destination + '/geom1.xyz'  
    shutil.copy(inpf_s_i, inpf_d_i)  
    shutil.copy(geom_s_i, geom_d_i)        
    if (to_add == 2):
        source = frontier_dirs[2]
        inpf_s_j = source + '/qommma.in'
        geom_s_j = source + '/geom1.xyz'
        destination = new_frontier_dirs[2]
        inpf_d_j = destination + '/qommma.in'
        geom_d_j = destination + '/geom1.xyz'  
        shutil.copy(inpf_s_j, inpf_d_j)  
        shutil.copy(geom_s_j, geom_d_j)   

    # Now, prepare the primitive internal coordinate displacements for qommma.in.
    prim_displacement_lst_i = []
    for indice,coord in enumerate(tangent_prims):
        dq = tangent[indice]
        to_append = (dq,) + (coord[0:4],)
        prim_displacement_lst_i.append(to_append)
    if (to_add==2):
        prim_displacement_lst_j = prim_displacement_lst_i.copy()
        for i in range(len(prim_displacement_lst_j)):
            prim_displacement_lst_j[i][0] = prim_displacement_lst_j[i][0] * (-1)

    # Now, prepare the new input prim_constrain_lst (which is the same if two new nodes are added) for qommma.in using the tangent and the primitive driving coordinates.
    prim_constrain_lst = []
    prim_constrain_lst_sub = []
    for indice,coord in enumerate(tangent_prims):
        to_append = coord[0:4]
        prim_constrain_lst_sub.append(to_append)
    prim_constrain_lst.append(prim_constrain_lst_sub)

    # Now, prepare the primitive internal coordinate linear combination coefficients.
    prim_coeff_lst_i = []
    prim_coeff_lst_sub = []
    for coeff in tangent:
        prim_coeff_lst_sub.append(coeff)  
    prim_coeff_lst_i.append(prim_coeff_lst_sub)
    if (to_add==2):
        prim_coeff_lst_j = prim_coeff_lst_i.copy()
        for i in range(len(prim_coeff_lst_j)):
            prim_coeff_lst_j[i] =  prim_coeff_lst_j[i] * (-1)
            
    
    # Now, the new qommma.in file in the new directories are modified to include the new tangent, the max number of steps, and the gsmphase.
    first_node_cons = True
    first_node_coeff = True
    first_node_disp = True
    first_node_disp_n = True
    first_node_type = True
    first_node_ncon = True
    first_node_phase = True
    first_node_maxcyc = True
    with open(inpf_d_i, 'r') as file:
        inp_data = file.readlines()
    for line_num,line in enumerate(inp_data):
        if "prim_constrain_lst" in line:
            inp_data[line_num] = "prim_constrain_lst=" + str(prim_constrain_lst) + '\n'
            first_node_cons = False
        if "prim_coeff_lst" in line:
            inp_data[line_num] = "prim_coeff_lst=" + str(prim_coeff_lst_i) + '\n'
            first_node_coeff = False
        if "prim_displacement_lst" in line:
            inp_data[line_num] = "prim_displacement_lst=" + str(prim_displacement_lst_i) + '\n'
            first_node_disp = False
        if "disp_prim" in line:
            inp_data[line_num] = "disp_prim=" + str(len(tangent)) + '\n'
            first_node_disp_n = False
        if "maxcycle" in line:
            inp_data[line_num] = 'maxcycle=5\n'
            first_node_maxcyc = False
        if "gsmphase" in line: 
            inp_data[line_num] = "gsmphase='growth'\n"
            first_node_phase = False
        if "gsmtype" in line:
            inp_data[line_num] = "gsmtype='se_gsm'\n"
            first_node_type = False
        if "ncon_prim" in line:
            inp_data[line_num] = "ncon_prim=" + str(1) + '\n'
            first_node_ncon = False
    with open(inpf_d, 'w') as file:
        file.writelines(inp_data)
        if first_node_coeff == True:
            file.write('\n')
            file.write("prim_coeff_lst=" + str(prim_coeff_lst_i) + '\n')
        if first_node_disp == True:
            file.write('\n')
            file.write("prim_displacement_lst=" + str(prim_displacement_lst_i) + '\n')
        if first_node_disp_n == True:
            file.write('\n')
            file.write("disp_prim=" + str(len(tangent)) + '\n')
        if first_node_cons == True:
            file.write('\n')
            file.write("prim_constrain_lst=" + str(prim_constrain_lst) + '\n')
        if first_node_type == True:
            file.write('\n')
            file.write("gsmtype='se_gsm'\n")
        if first_node_maxcyc == True:
            file.write('\n')
            file.write('maxcycle=5\n')
        if first_node_ncon == True:
            file.write('\n')
            file.write("ncon_prim=" + str(1) + '\n')            
        if first_node_phase == True:
            file.write('\n')
            file.write("gsmphase='growth'\n")
    if (to_add==2):
        first_node_cons = True
        first_node_coeff = True
        first_node_disp = True
        first_node_disp_n = True
        first_node_type = True
        first_node_ncon = True
        first_node_phase = True
        first_node_maxcyc = True
        with open(inpf_d_j, 'r') as file:
            inp_data = file.readlines()
        for line_num,line in enumerate(inp_data):
            if "prim_constrain_lst" in line:
                inp_data[line_num] = "prim_constrain_lst=" + str(prim_constrain_lst) + '\n'
                first_node_cons = False
            if "prim_coeff_lst" in line:
                inp_data[line_num] = "prim_coeff_lst=" + str(prim_coeff_lst_j) + '\n'
                first_node_coeff = False
            if "prim_displacement_lst" in line:
                inp_data[line_num] = "prim_displacement_lst=" + str(prim_displacement_lst_j) + '\n'
                first_node_disp = False
            if "disp_prim" in line:
                inp_data[line_num] = "disp_prim=" + str(len(tangent)) + '\n'
                first_node_disp_n = False
            if "maxcycle" in line:
                inp_data[line_num] = 'maxcycle=5\n'
                first_node_maxcyc = False
            if "gsmphase" in line: 
                inp_data[line_num] = "gsmphase='growth'\n"
                first_node_phase = False
            if "gsmtype" in line:
                inp_data[line_num] = "gsmtype='de_gsm'\n"
                first_node_type = False
            if "ncon_prim" in line:
                inp_data[line_num] = "ncon_prim=" + str(1) + '\n'
                first_node_ncon = False
        with open(inpf_d, 'w') as file:
            file.writelines(inp_data)
            if first_node_coeff == True:
                file.write('\n')
                file.write("prim_coeff_lst=" + str(prim_coeff_lst_j) + '\n')
            if first_node_disp == True:
                file.write('\n')
                file.write("prim_displacement_lst=" + str(prim_displacement_lst_j) + '\n')
            if first_node_disp_n == True:
                file.write('\n')
                file.write("disp_prim=" + str(len(tangent)) + '\n')
            if first_node_cons == True:
                file.write('\n')
                file.write("prim_constrain_lst=" + str(prim_constrain_lst) + '\n')
            if first_node_type == True:
                file.write('\n')
                file.write("gsmtype='de_gsm'\n")
            if first_node_maxcyc == True:
                file.write('\n')
                file.write('maxcycle=5\n')
            if first_node_ncon == True:
                file.write('\n')
                file.write("ncon_prim=" + str(1) + '\n')            
            if first_node_phase == True:
                file.write('\n')
                file.write("gsmphase='growth'\n")   
            
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
    
    # NOTE: need to code this...
    #returns nothing, but makes new files for the reparameterisation and calls fortran code to update the geometry.

def DE_initialise_prims(prim_list_R, prim_list_P, usrdir):
    """
    
    // Function which reads the set of primitive internal coordinates from the nodeR and nodeP directories which describe the reaction pathway. //
    // This function is used in the set-up of the double-ended GSM. //
    
    Arguments
    ----------
    prim_list_R : string
        Directory to the file prim_list in the reactant node directory.
    prim_list_P : string
        Directory to the file prim_list in the product node directory.
    usrdir : string
        The user directory.

    """
    
    # NOTE: need to fix so that it reads both reactant and product primitives, and then makes sure that they are the same.
    # If they are not the same, then end the program and suggest either running nodeR and nodeP with xyz additional primitives.
    # Or, suggest that they use total connectivity instead.

    # Opening the file...
    with open(prim_list_R, 'r') as r_f:
        # First, read the number of primitive internal coordinates from the top of the file.
        prim_num = int(r_f.readline())
        
        # Now, read all the primitive internal coordinate definitions.
        prim_defs = []
        for i in range(0,prim_num):
            temp_lst = []
            temp_str = r_f.readline()
            temp_lst = temp_str.split()
            prim_defs.append(temp_lst)
    
    return prim_num, prim_defs   

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
    for indice,drive in enumerate(driving_coords):
        driving = drive[0:4]
        direction = drive[4]

        # Now, find out what type of primitive internal coordinate each driving coordinate is by counting the number of zeroes.
        zero_count = 0
        for j in driving:
            if (j == 0):
                zero_count += 1
        
        if (zero_count == 2): # Bond
            tangent.append(0.1 * direction) # Angstroms...
        elif (zero_count == 1): # Angle
            tangent.append(0.0872665 * direction) # Radians (5 deg)...
        elif (zero_count == 0): # Dihedral torsion
            tangent.append(0.0872665 * direction) # Radians (5 deg)...
           
    return tangent

def SE_get_final_tangent(current_nodes, driving_coords, usrdir):
    """
    
    // Function which generates the primitive internal coordinate driving coordinate which describes the reaction pathway. //
    // This function is only used during the growth-phase of SE-GSM, and specifically for the final 2 - 4 nodes. //
    // The magnitude of the tangent is scaled based on how much of a difference is found between the frontier node and the reactant node. //
    // This is approximate, and assumes that the energy profile is roughly symmetrical (not usually true), but works for growing the string. //
    
    Arguments
    ----------
    current_nodes : integer
        The current number of nodes, including the reactant directory.
    driving_coords : list
        The primitive internal coordinates which the user has specified to be likely involved in the reaction pathway.
    usrdir : string
        The user directory.

    """
    
    # This is relatively similar to the function SE_get_tangent, except that it returns a total tangent for 2-4 nodes.
    tangent = []
    for indice,drive in enumerate(driving_coords):
        driving = drive[0:4]
        direction = drive[4]

        # Now, find out what type of primitive internal coordinate each driving coordinate is by counting the number of zeroes.
        zero_count = 0
        for j in driving:
            if (j == 0):
                zero_count += 1
        
        if (zero_count == 2): # Bond
            tangent.append(0.1 * direction)# * current_nodes) # Angstroms...
        elif (zero_count == 1): # Angle
            tangent.append(0.0872665 * direction)# * current_nodes) # Radians (5 deg)...
        elif (zero_count == 0): # Dihedral torsion
            tangent.append(0.0872665 * direction)# * current_nodes) # Radians (5 deg)...
    
    # Make the reaction pathway symmetrical.
    total_new = current_nodes
    return tangent, total_new
              
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
    geom_s = source + '/geom1.xyz'
    destination = new_frontier_dir
    inpf_d = destination + '/qommma.in'
    geom_d = destination + '/init_geom1.xyz'  
    shutil.copy(inpf_s, inpf_d)  
    shutil.copy(geom_s, geom_d)  
    
    # Now, prepare the primitive internal coordinate displacements for qommma.in.
    prim_displacement_lst = []
    for indice,drive in enumerate(driving_coords):
        dq = tangent[indice]
        to_append = (dq,) + (drive[0:4],)
        prim_displacement_lst.append(to_append)
        
    # Now, prepare the new input prim_constrain_lst for qommma.in using the tangent and the primitive driving coordinates.
    prim_constrain_lst = []
    prim_constrain_lst_sub = []
    for indice,drive in enumerate(driving_coords):
        to_append = drive[0:4]
        prim_constrain_lst_sub.append(to_append)
    prim_constrain_lst.append(prim_constrain_lst_sub)
    
    # Now, prepare the primitive internal coordinate linear combination coefficients. This is simple for the growth phase.
    prim_coeff_lst = []
    prim_coeff_lst_sub = []
    for indice,dq in enumerate(tangent):
        prim_coeff_lst_sub.append(math.copysign(1,dq))
    prim_coeff_lst.append(prim_coeff_lst_sub)
    
    # Now, the new qommma.in file in the new directory is modified to include the new tangent, the max number of steps, and the gsmphase.
    first_node_cons = True
    first_node_coeff = True
    first_node_disp = True
    first_node_disp_n = True
    first_node_type = True
    first_node_ncon = True
    first_node_phase = True
    first_node_maxcyc = True
    with open(inpf_d, 'r') as file:
        inp_data = file.readlines()
    for line_num,line in enumerate(inp_data):
        if "prim_constrain_lst" in line:
            inp_data[line_num] = "prim_constrain_lst=" + str(prim_constrain_lst) + '\n'
            first_node_cons = False
        if "prim_coeff_lst" in line:
            inp_data[line_num] = "prim_coeff_lst=" + str(prim_coeff_lst) + '\n'
            first_node_coeff = False
        if "prim_displacement_lst" in line:
            inp_data[line_num] = "prim_displacement_lst=" + str(prim_displacement_lst) + '\n'
            first_node_disp = False
        if "disp_prim" in line:
            inp_data[line_num] = "disp_prim=" + str(len(tangent)) + '\n'
            first_node_disp_n = False
        if "maxcycle" in line:
            inp_data[line_num] = 'maxcycle=5\n'
            first_node_maxcyc = False
        if "gsmphase" in line: 
            inp_data[line_num] = "gsmphase='growth'\n"
            first_node_phase = False
        if "gsmtype" in line:
            inp_data[line_num] = "gsmtype='se_gsm'\n"
            first_node_type = False
        if "ncon_prim" in line:
            inp_data[line_num] = "ncon_prim=" + str(1) + '\n'
            first_node_ncon = False
    with open(inpf_d, 'w') as file:
        file.writelines(inp_data)
        if first_node_coeff == True:
            file.write('\n')
            file.write("prim_coeff_lst=" + str(prim_coeff_lst) + '\n')
        if first_node_disp == True:
            file.write('\n')
            file.write("prim_displacement_lst=" + str(prim_displacement_lst) + '\n')
        if first_node_disp_n == True:
            file.write('\n')
            file.write("disp_prim=" + str(len(tangent)) + '\n')
        if first_node_cons == True:
            file.write('\n')
            file.write("prim_constrain_lst=" + str(prim_constrain_lst) + '\n')
        if first_node_type == True:
            file.write('\n')
            file.write("gsmtype='se_gsm'\n")
        if first_node_maxcyc == True:
            file.write('\n')
            file.write('maxcycle=5\n')
        if first_node_ncon == True:
            file.write('\n')
            file.write("ncon_prim=" + str(1) + '\n')            
        if first_node_phase == True:
            file.write('\n')
            file.write("gsmphase='growth'\n")        
       
def SE_add_final_nodes(frontier_dir, new_frontier_dirs, tangent, driving_coords, usrdir):
    """
    
    // Function which creates new qommma.in files for the final nodes which essentially represent the new frontier nodes. //
    // The geometry from the previous frontier node is used as the starting point for all the new nodes and also copied over to the new directories. //
    
    Arguments
    ----------
    frontier_dir : string
        String of the current frontier node which is about to be replaced.
    new_frontier_dirs : list
        List of the 2 - 4 new frontier node directories which are about to be added.
    tangent : list
        The primitive internal coordinate tangent which is used to define the geometry and constraints applied to the new nodes.
    driving_coords : list
        The primitive internal coordinates which the user has specified to be likely involved in the reaction pathway.
    usrdir : string
        The user directory.

    """
    
    # First, initialise the source directory.
    source = frontier_dir
    inpf_s = source + '/qommma.in'
    
    # Now, iterate through each of the new frontier node directories and copy the files appropriately.
    for dir in new_frontier_dirs:
        # First, prepare the directory and copy over files.
        destination = dir
        inpf_d = destination + '/qommma.in'
        shutil.copy(inpf_s, inpf_d)   
            
        # Now, prepare the primitive internal coordinate displacements for qommma.in.
        prim_displacement_lst = []
        for indice,drive in enumerate(driving_coords):
            dq = tangent[indice]
            to_append = (dq,) + (drive[0:4],)
            prim_displacement_lst.append(to_append)
            
        # Now, prepare the new input prim_constrain_lst for qommma.in using the tangent and the primitive driving coordinates.
        prim_constrain_lst = []
        prim_constrain_lst_sub = []
        for indice,drive in enumerate(driving_coords):
            to_append = drive[0:4]
            prim_constrain_lst_sub.append(to_append)
        prim_constrain_lst.append(prim_constrain_lst_sub)
        
        # Now, prepare the primitive internal coordinate linear combination coefficients. This is simple for the growth phase.
        prim_coeff_lst = []
        prim_coeff_lst_sub = []
        for dq in tangent:
            prim_coeff_lst_sub.append(math.copysign(1,dq))
        prim_coeff_lst.append(prim_coeff_lst_sub)
          
        # Now, the new qommma.in file in the new directory is modified to include the new tangent, the max number of steps, and the gsmphase.
        with open(inpf_d, 'r') as file:
            inp_data = file.readlines()
        for line_num,line in enumerate(inp_data):
            if "prim_constrain_lst" in line:
                inp_data[line_num] = "prim_constrain_lst=" + str(prim_constrain_lst) + '\n'
            if "prim_coeff_lst" in line:
                inp_data[line_num] = "prim_coeff_lst=" + str(prim_coeff_lst) + '\n'
            if "prim_displacement_lst" in line:
                inp_data[line_num] = "prim_displacement_lst=" + str(prim_displacement_lst) + '\n'
            if "disp_prim" in line:
                inp_data[line_num] = "disp_prim=" + str(len(tangent)) + '\n'
            if "maxcycle" in line:
                inp_data[line_num] = 'maxcycle=5\n'
            if "gsmphase" in line: 
                inp_data[line_num] = "gsmphase='growth'\n"
            if "gsmtype" in line:
                inp_data[line_num] = "gsmtype='se_gsm'\n"
            if "ncon_prim" in line:
                inp_data[line_num] = "ncon_prim=" + str(1) + '\n'
        with open(inpf_d, 'w') as file:
            file.writelines(inp_data)
    
def SE_check_delE(old_frontier_dir, frontier_dir, current_nodes):
    """
    
    // Function which checks the difference in energy between the current and previous frontier nodes. //
    // If the energy of the new frontier node is lower, then a stationary point has been passed and we can move on. //
    
    Arguments
    ----------
    old_frontier_dir : string
        String of the old frontier node directory.
    frontier_dir : string
        String of the current frontier node directory.
    current_nodes : integer
        The current number of nodes, including the reactant directory.

    """
    
    # The relevant energies are most easily found in the report files.
    old_report = old_frontier_dir + '/report1'
    cur_report = frontier_dir + '/report1'

    # Now, open the reports and reverse the order to find the last energy evaluation.
    with open(old_report, 'r') as old_f:
        lines = old_f.readlines()
        for line in reversed(lines):
            if 'Present, previous energy' in line:
                split_line = line.split()
                old_E = float(split_line[3])
                break
    with open(cur_report, 'r') as cur_f:
        lines = cur_f.readlines()
        for line in reversed(lines):
            if 'Present, previous energy' in line:
                split_line = line.split()
                cur_E = float(split_line[3])
                break
    
    # Now, evaluate the difference and decide if a stationary point has been surpassed.
    del_E = cur_E - old_E
    if ((del_E < 0.0) and (current_nodes > 8)):
        SP_found = True
    else:
        SP_found = False

    return SP_found
    
def SE_initialise_prims(prim_list_R, usrdir):
    """
    
    // Function which reads the set of primitive internal coordinates from the nodeR directory which describe the reaction pathway. //
    // This function is used in the set-up of the single-ended GSM. //
    
    Arguments
    ----------
    prim_list_R : string
        Directory to the file prim_list in the reactant node directory.
    usrdir : string
        The user directory.

    """
    
    # Opening the file...
    with open(prim_list_R, 'r') as r_f:
        # First, read the number of primitive internal coordinates from the top of the file.
        prim_num = int(r_f.readline())
        
        # Now, read all the primitive internal coordinate definitions.
        prim_defs = []
        for i in range(0,prim_num):
            temp_lst = []
            temp_str = r_f.readline()
            temp_lst = temp_str.split()
            prim_defs.append(temp_lst)
    
    return prim_num, prim_defs

def additional_prims(add_prims, prim_num, prim_defs, usrdir):
    """
    
    // Function which adds any additional primitive internal coordinate definitions to the global definitions. //
    
    Arguments
    ----------
    add_prims : list
        List containing the additional primitive internal coordinate definitions to be added.
    prim_num : integer
        The total number of primitive internal coordinates.
    prim_defs : list
        List containing the primitive internal coordinate definitions.
    usrdir : string
        The user directory.

    """
    
    # Simply iterate through the additional primitives, add them to the primitive definitions, and update the number of primitives.
    new_prims = 0
    for i in add_prims:
        prim_defs.append(i)
        new_prims += 1
    prim_num += new_prims
           
    return prim_num, prim_defs
    
def extract_energies(node_dirs):
    """
    
    // Function which extracts the energies of all nodes and adds to a list. //
    
    Arguments
    ----------
    node_dirs : list
        List of all node directories.
    frontier_dir : string
        String of the current frontier node directory.

    """
    
    # Iterate through all the node directories.
    energies = []
    for dir in node_dirs:
        # The relevant energies are most easily found in the report files.
        report = dir + '/report1'

        # Now, open the report in reverse order to find the last energy evaluation.
        with open(report, 'r') as f:
            lines = f.readlines()
            for line in reversed(lines):
                if 'Present, previous energy' in line:
                    split_line = line.split()
                    E = float(split_line[3])
                    break
        
        # Lastly, add to the list...
        energies.append(E)

    return energies
  
def get_tangents_opt(node_dirs, usrdir, driving_coords = None):
    """
    
    // Function which generates the set of primitive internal coordinate tangents which describe the reaction pathway. //
    // This function is only used during the optimisation-phase of GSM, and for both the single- and double-ended variant. //
    
    Arguments
    ----------
    node_dirs : list
        List of all node directories which are used to define the tangents.
    usrdir : string
        The user directory.
    driving_coords : list
        The primitive internal coordinates which the user has specified to be likely involved in the reaction pathway.
        Implemented as a kwarg for the case of the double-ended method.

    """
    
    # These tangents are based on Henkelman's technique to minimise kinks in the string.
    # Simply, tangents are defined to be pointing in the direction of highest energy.
    # Thus, we need to first get the energies of all nodes.
    energies = extract_energies(node_dirs)

    # Iterate through each of the directories.
    tangent_list = []
    tangent_prims_list = [] 
    for counter,dir_i in enumerate(node_dirs):
        if not ("nodeR" or "nodeP") in dir_i:
            # First, read in the primitive internal coordinates from the ith node directory.
            dir_ii = dir_i + '/jobfiles/'
            prims_i = []
            prim_list_i = []
            os.chdir(dir_ii)
            with open("prims", 'r') as pr, open("prim_list", 'r') as pr_ls:
                n_prims = int(pr.readline())
                pr_ls.readline()
                for i in range(0,n_prims):
                    prims_i.append(float(pr.readline()))
                    prim_list_i.append(pr_ls.readline())
            
            # Now, check the energy of the node and compare to neighbouring nodes.
            # Also account for the special cases of the first and last node.
            if counter == 0:
                dir_j = node_dirs[1]
            elif counter == (len(node_dirs) - 1):
                dir_j = node_dirs[len(node_dirs)-2]
            else:
                energy_imin1 = energies[counter-1]
                energy_iplus1 = energies[counter+1]
                if (energy_imin1 > energy_iplus1):
                    dir_j = node_dirs[counter-1]
                else:
                    dir_j = node_dirs[counter+1]
            
            # Read in the primitive internal coordinates from the jth node directory.
            dir_jj = dir_j + '/jobfiles/'
            prims_j = []
            prim_list_j = []
            os.chdir(dir_jj)
            with open("prims", 'r') as pr, open("prim_list", 'r') as pr_ls:
                n_prims = int(pr.readline())
                pr_ls.readline()
                for i in range(0,n_prims):
                    prims_j.append(float(pr.readline()))
                    prim_list_j.append(pr_ls.readline()) 

            # Now, the tangent can be calculated simply as the difference between these primitive internal coordinates.
        tangent_temp = []
        for p1, p2 in zip(prims_i, prims_j):
            p = p1 - p2
            tangent_temp.append(p)
            
            # If a given primitive is below a certain threshold, then it is set to zero and hence not constrained in the Fortran optimisation.
            for i in range(0,n_prims):
                # Store the tangent and primitive associated with it.
                temp_tan = tangent_temp[i]
                temp_prim = prim_list_i[i].split()
                
                # Store the type of primitive internal coordinate.
                zero_count = 0
                for prim in temp_prim:
                    if int(prim) == 0:
                        zero_count += 1
                
                # Setting appropriate element to zero.
                if (abs(temp_tan) < 1E-2):
                    tangent_temp[i] = 0.0  
   
                ################################################################################
                # TEMPORARY TEST FIX - means that valence and dihedral angles are not included.#
                ################################################################################
                if (zero_count == 0) or (zero_count == 1):
                    tangent_temp[i] = 0.0

            # We should remove components of the tangent which are zero and save the relevant primitive definitions corresponding to the tangents.
            tangent = []
            tangent_prims = []
            for counter,tang in enumerate(tangent_temp):
                if abs(tang) > 0.0:
                    tangent.append(tang)
                    tangent_prims.append(prim_list_i[counter])      
                    
            # Lastly, append both the tangent and the list of relevant primitives to tangent_list and tangent_prims_list.
            tangent_list.append(tangent)
            tangent_prims_list.append(tangent_prims)
        else:
            tangent = []
            tangent_prims = []
            tangent_list.append(tangent)
            tangent_prims_list.append(tangent_prims)

    return tangent_list, tangent_prims_list
    
def gen_input_opt(node_dirs, tangent_list, tangent_prims_list, usrdir):
    """
    
    // Function which creates a series of new qommma.in files for all the nodes. //
    // This function is only used during the optimisation-phase of GSM, and for both the single- and double-ended variant. //
    
    Arguments
    ----------
    node_dirs : list
        List of all node directories which are used to define the tangents.
    tangent_list : array
        Array of all the primitive internal coordinate tangents which connect the nodes and create the reaction pathway.
    tangent_prims_list : array
        Array of all the primitive internal coordinate which define the tangents.
    usrdir : string
        The user directory.

    """

    # Now, iterate through each of the new frontier node directories and copy the files appropriately.
    for counter,dir in enumerate(node_dirs):
        if not ("nodeR" or "nodeP") in dir:
            # First, initialise the source and destination directories (which are the same).
            source = dir
            destination = dir
            
            # Now copy the geometry from the most recent QoMMMa optimisation.
            geom_s = source + '/geom1.xyz'
            geom_d = destination + '/init_geom1.xyz'
            shutil.copy(geom_s, geom_d)
            os.remove(geom_s)

            # Prepare the new input prim_constrain_lst for the given qommma.in using the tangent and the primitive internal coordinates.
            prim_constrain_lst = []
            prim_constrain_lst_sub = []
            for prim in tangent_prims_list[counter]:
                prim.replace("\n","")
                prim = prim.split()
                to_append = prim
                prim_constrain_lst_sub.append(to_append)
            prim_constrain_lst.append(prim_constrain_lst_sub)
            
            # Now, prepare the primitive internal coordinate linear combination coefficients.
            prim_coeff_lst = []
            prim_coeff_lst_sub = []
            for coeff in tangent_list[counter]:
                prim_coeff_lst_sub.append(coeff)  
            prim_coeff_lst.append(prim_coeff_lst_sub)
            
            # Now, the new qommma.in file in the new directory is modified to include the new tangent, the max number of steps, and the gsmphase.
            inpf_d = destination + '/qommma.in'
            with open(inpf_d, 'r') as file:
                inp_data = file.readlines()
            for line_num,line in enumerate(inp_data):
                if "prim_constrain_lst" in line:
                    inp_data[line_num] = "prim_constrain_lst=" + str(prim_constrain_lst) + '\n'
                if "prim_coeff_lst" in line:
                    inp_data[line_num] = "prim_coeff_lst=" + str(prim_coeff_lst) + '\n'
                if "prim_displacement_lst" in line:
                    inp_data[line_num] = "prim_displacement_lst='None'\n"
                if "disp_prim" in line:
                    inp_data[line_num] = "disp_prim=0\n"
                if "maxcycle" in line:
                    inp_data[line_num] = 'maxcycle=15\n'
                if "gsmphase" in line: 
                    inp_data[line_num] = "gsmphase='opt'\n"
                if "gsmtype" in line:
                    inp_data[line_num] = "gsmtype='se_gsm'\n"
                if "ncon_prim" in line:
                    inp_data[line_num] = "ncon_prim=" + str(1) + '\n'
            with open(inpf_d, 'w') as file:
                file.writelines(inp_data)

def check_convergence(node_dirs, current_nodes, usrdir):
    """
    
    // Function which checks the convergence of all nodes. //
    // If all nodes have converged, then the function returns True. //
    
    Arguments
    ----------
    node_dirs : list
        List of all node directories on the string.
    usrdir : string
        The user directory.

    """
    
    # Initialising the convergence boolean logic array.
    is_converged = [False] * current_nodes
    
    for counter,dir in enumerate(node_dirs):
        # The convergence is most easily found in the report files.
        report = dir + '/report1'

        # Now, open the reports and reverse the order to find the relevant string.
        with open(report, 'r') as f:
            lines = f.readlines()
            for line in reversed(lines):
                if 'Congratulations!!' in line:
                    is_converged[counter] = True
                    break

    return is_converged
    
def reparam_opt(node_dirs, total_nodes, usrdir):
    """
    
    // Function which reparameterises the string. //
    // This is to say that the geometries of all nodes are evenly spaced along the tangent defined between reactant and product node. //
    
    Arguments
    ----------
    node_dirs : list
        List of all node directories on the string.
    total_nodes : integer
        Integer representing the total number of nodes on the string.
    usrdir : string
        The user directory.

    """
    
    # To equally respace the nodes, the tangent between reactant and product nodes must first be obtained.
    nodeR_dir = node_dirs[0]
    nodeP_dir = node_dirs[-1]

    # Read the primitive definitions and values for the reactant node.
    dir_r = nodeR_dir + '/jobfiles/'
    prims_r = []
    prim_list_r = []
    os.chdir(dir_r)
    with open("prims", 'r') as pr, open("prim_list", 'r') as pr_ls:
        n_prims = int(pr.readline())
        pr_ls.readline()
        for i in range(0,n_prims):
            prims_r.append(float(pr.readline()))
            prim_list_r.append(pr_ls.readline())  
    
    # Now for the product node.
    dir_p = nodeP_dir + '/jobfiles/'
    prims_p = []
    prim_list_p = []
    os.chdir(dir_p)
    with open("prims", 'r') as pr, open("prim_list", 'r') as pr_ls:
        n_prims = int(pr.readline())
        pr_ls.readline()
        for i in range(0,n_prims):
            prims_p.append(float(pr.readline()))
            prim_list_p.append(pr_ls.readline())  

    # Find the overall tangent which the nodes will be spaced along.
    # If the change in primitive internal coordinates is not particularly great, then do not space the node along this coordinate.
    total_tangent = []
    total_tangent_prims = []
    for i in range(0,n_prims):
        p1 = prims_r[i]
        p2 = prims_p[i]
        p = p2 - p1

        # Store the type of primitive internal coordinate.
        temp_prim = prim_list_r[i].split()
        zero_count = 0
        for prim in temp_prim:
            if int(prim) == 0:
                zero_count += 1
        
        # Populate total_tangent appropriately.
        if (zero_count == 2):
            total_tangent.append(p)  
        else:
            total_tangent.append(0.0)
        total_tangent_prims.append(prim_list_r[i])
 
    # Initialise the starting point for the primitive internal coordinates which the reparameterisation will occur across.
    # This array is maintained as the generation of qommma.in files proceeds.
    starting_prims = prims_r
        
    # Now, new qommma.in files are created for each of the node directories so that the string can be reparameterised.
    for counter,dir in enumerate(node_dirs):
        if not ("nodeR" or "nodeP") in dir:
            # First, initialise the source and destination directories (which are the same).
            source = dir
            destination = dir
            
            # Now, get the primitive internal coordinates of the ith node...
            dir_i = dir + '/jobfiles/'
            prims_i = []
            os.chdir(dir_i)
            with open("prims", 'r') as pr:
                    n_prims = int(pr.readline())
                    for i in range(0,n_prims):
                        prims_i.append(float(pr.readline()))
           
            # Generate the primitive internal coordinates from reparameterisation.
            reparam_prims = []
            for i in range(0,n_prims):
                reparam_prims.append(starting_prims[i] + (total_tangent[i] / total_nodes))
                
            # Update starting prims...
            starting_prims = reparam_prims
            
            # Now, compare the actual primitives (prims_i) and the primitives generated by reparameterisation (reparam_prims).
            # If the difference between them is large enough, then include in the reparameterisation.
            tangent_prims_list = []
            tangent_list = []
            for i in range(0,n_prims):
                p1 = prims_i[i]
                p2 = reparam_prims[i]
                p = p2 - p1
                if (abs(p) > 0.01) and (abs(total_tangent[i]) > 0.0):
                    tangent_prims_list.append(prim_list_r[i])
                    tangent_list.append(p)    
     
            # Prepare the new input prim_constrain_lst for the given qommma.in using the tangent and the primitive internal coordinates.
            # In this case, there is no step to be taken because we are in the optimisation phase.
            prim_displacement_lst = []
            for prim,tang in zip(tangent_prims_list,tangent_list):
                prim.replace("\n","")
                prim = prim.split()
                dq = tang
                to_append = (dq,) + (prim,)
                prim_displacement_lst.append(to_append)
            
            # Now, the new qommma.in file in the new directory is modified to include the new tangent and the gsmphase.
            inpf_d = destination + '/qommma.in'
            with open(inpf_d, 'r') as file:
                inp_data = file.readlines()
            for line_num,line in enumerate(inp_data):
                if "prim_constrain_lst" in line:
                    inp_data[line_num] = "prim_constrain_lst='None'" + '\n'
                if "prim_coeff_lst" in line:
                    inp_data[line_num] = "prim_coeff_lst='None'" + '\n'
                if "prim_displacement_lst" in line:
                    inp_data[line_num] = "prim_displacement_lst=" + str(prim_displacement_lst) + '\n'
                if "disp_prim" in line:
                    inp_data[line_num] = "disp_prim=" + str(len(tangent_list)) + '\n'
                if "ncon_prim" in line:
                    inp_data[line_num] = "ncon_prim=" + str(0) + '\n'
            with open(inpf_d, 'w') as file:
                file.writelines(inp_data) 
        

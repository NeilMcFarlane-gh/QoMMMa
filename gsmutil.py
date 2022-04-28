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
    
def gsmend(s, cwd, usrdir):
    """
    
    // Function which ends the given GSM QoMMMa calculation and writes the error to the log file QoMMMa8_GSM.log. //
    
    Arguments
    ----------
    s : string
        The message to be written to the file QoMMMa8.log.
    cwd : string
        The current working directory.
    usrdir : string
        The user directory.

    """
 
    # The time the program ends is recorded along with an error message.
    hlp = asctime()
    qomlog(s + '\n' + '  QoMMMa job ends at ' + str(hlp), usrdir)
    
    # Ending QoMMMa GSM job.
    sys.exit()


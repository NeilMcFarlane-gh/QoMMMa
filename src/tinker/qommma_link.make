#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##                  (Linux/GNU g77 Version)                  ##
#  ##                                                           ##
#  ###############################################################
#
	gfortran-6  -s -o analyze_qommma.x analyze_qommma.o libtinker.a
	gfortran-6  -s -o minimize_qommma.x minimize.o libtinker.a

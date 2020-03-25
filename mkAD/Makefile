ifeq "${TOOL}" "tap"
    LOG_MSG=
    copy_extra=cp solve_dummy.f90 ${PROCESS_DIR}${dir}/solve.f90
endif

FCALL_STR=f${TOOL}_sub
RCALL_STR=r${TOOL}_sub
SHEMAT_DIR=..
FINAL_DESTINY=${SHEMAT_DIR}/g_${TOOL}/${phys_base}
FINAL_DESTINY_ADJ=${SHEMAT_DIR}/ad_${TOOL}/${phys_base}
PROCESS_DIR=process_${TOOL}
files_filter=*.f *.f77 *.f90 *.f95 *.F *.inc
PP_FLAGS_base =-DG95 -DCLopt -DUSE_QDDOT
#PP_FLAGS ="-DG95 -DCLopt -DQDDOT $(addprefix -DPROPS_,${props_lst}) $(addprefix -DUSER_,${users_lst})"
ifeq "${phys_base}" "head"
    PP_FLAGS_base += -Dhead_base
endif
ifeq "${phys_base}" "pres"
    PP_FLAGS_base += -Dpres_base
endif

#
# !!! for testing only !!! normaly defined outside !!!
#props_lst=basc
#users_lst=wells3d
#props_lst=const
#user_lst=none
#forward_lst=hdf5 solve forward/head forward blas box nonlinear
# !!! -------------------------------------------- !!!

export AD_MAKE_DIR=${CURDIR}
DESTINY_alts=$(foreach udir, ${users_lst}, $(foreach pdir, ${props_lst}, /${pdir}/${udir}))
DESTINY_alts_mkdir= $(foreach pdir, ${props_lst}, /${pdir})

SHELL:=/bin/bash

default:
	@echo 'type:                     $(MAKE)   <target>'
	@echo 'where target is one of:     prepare | tlm | adm | clean '
	@echo 'with the following meanings: '
	@echo ' prepare = Copy and preprocess files and remove unneeded ones in all subdirectories.'
	@echo ' prepare_reverse = Copy and preprocess files and remove unneeded ones in all subdirectories (reverse mode).'
	@echo ' tlm     = Run AD-Tool (Tapenade) to generate augmented FORTRAN-codes for derivatives in all subdirectories (tangent linear = gradient).'
	@echo ' adm     = Run AD-Tool (Tapenade) to generate augmented FORTRAN-codes for derivatives  in all subdirectories (adjoint`¸s).'
	@echo ' clean   = objects, AD-generated code, and executables.'

setup:
	@echo Copying files from the following directories:
	@echo "  " ${forward_lst}
	@echo to $(words ${DESTINY_alts}) subdirectories [$(DESTINY_alts)] for the alternative combinations
	#clear the previously created folders and recreate them
	@rm -rf ${PROCESS_DIR}
	@mkdir ${PROCESS_DIR}
	-@$(foreach dir, ${DESTINY_alts}, rm -rf ${PROCESS_DIR}${dir};) rm -rf ${PROCESS_DIR}/tmp;
	#create new folders based on DESTINY_alts e.g. ${PROCESS_DIR}/const/none and ${PROCESS_DIR}/tmp
	@echo "Creating subdirectories in ${PROCESS_DIR}..."; $(foreach dir, ${DESTINY_alts}, echo "creating ${PROCESS_DIR}${dir}";mkdir -p ${PROCESS_DIR}${dir};) mkdir ${PROCESS_DIR}/tmp
	#copy all forward-files to be transformed into tmp
	@$(foreach dir, ${forward_lst}, cp $(wildcard $(addprefix ${SHEMAT_DIR}/${dir}/,${files_filter})) ${PROCESS_DIR}/tmp/;)
	#copy the Makefile.sub into tmp
	@cp Makefile.sub ${PROCESS_DIR}/tmp/Makefile
	#copy everything from tmp into the ${DESTINY_alts} directories
	@$(foreach dir, ${DESTINY_alts}, cp ${PROCESS_DIR}/tmp/*  ${PROCESS_DIR}${dir}/;)
	#copy all user/props files into the appropriate folders (e.g. fortran-files from user/none and props/const into ${PROCESS_DIR}/const/none)
	@$(foreach diru, ${users_lst}, $(foreach dirp, ${props_lst}, cp $(wildcard $(addprefix ${SHEMAT_DIR}/user/${diru}/,${files_filter}) $(addprefix ${SHEMAT_DIR}/props/${dirp}/,${files_filter})) ${PROCESS_DIR}/${dirp}/${diru}/;))

prepare: setup
#	@$(MAKE) setup TOOL=${TOOL} phys_base=${phys_base} ad_mode=${ad_mode}
	@$(foreach dir, ${DESTINY_alts}, cd ${PROCESS_DIR}${dir};$(MAKE) setup_cleanup phys_base=${phys_base} ad_mode=${ad_mode}; cd ${CURDIR};)
	@$(foreach dir, ${DESTINY_alts}, ${copy_extra}; cd ${PROCESS_DIR}${dir}; $(MAKE) setup_dummycalls TOOL=${TOOL} phys_base=${phys_base} ad_mode=${ad_mode}; cd ${CURDIR};)

prepare_reverse: setup
#	@$(MAKE) setup TOOL=${TOOL} phys_base=${phys_base} ad_mode=${ad_mode}
	@$(foreach dir, ${DESTINY_alts}, cd ${PROCESS_DIR}${dir};$(MAKE) setup_cleanup phys_base=${phys_base} ad_mode=${ad_mode}; cd ${CURDIR};)
	@$(foreach dir, ${DESTINY_alts}, ${copy_extra}; cd ${PROCESS_DIR}${dir}; $(MAKE) setup_dummycalls_reverse TOOL=${TOOL} phys_base=${phys_base} ad_mode=${ad_mode} PP_FLAGS="${PP_FLAGS_base} -DPROPS_$(subst /,,$(dir ${dir})) -DUSER_$(notdir ${dir})"; cd ${CURDIR};)

tlm:
	@$(foreach dir, ${DESTINY_alts}, cd ${PROCESS_DIR}${dir}; $(MAKE) ${FCALL_STR} TOOL=${TOOL} phys_base=${phys_base} ad_mode=${ad_mode}; cd ${CURDIR};)
	-@$(foreach dir, ${DESTINY_alts}, rm -rf ${FINAL_DESTINY}${dir};)
	@echo "Creating subdirectories in ${FINAL_DESTINY}..."; $(foreach dir, ${DESTINY_alts}, mkdir -p ${FINAL_DESTINY}${dir};)
	-@$(foreach dir, ${DESTINY_alts}, cp ${PROCESS_DIR}${dir}/grad/* ${LOG_MSG} ${FINAL_DESTINY}${dir}/;)
	-@mkdir -p ${AD_MAKE_DIR}/../user/none/g_${TOOL}/${phys_base}
	@cp ${AD_MAKE_DIR}/g_calc_user.f90 ${AD_MAKE_DIR}/../user/none/g_${TOOL}/${phys_base}/.
	@$(MAKE) showall TOOL=${TOOL} phys_base=${phys_base} ad_mode=${ad_mode}

adm:
	@$(foreach dir, ${DESTINY_alts}, cd ${PROCESS_DIR}${dir}; $(MAKE) ${RCALL_STR} TOOL=${TOOL} phys_base=${phys_base} ad_mode=${ad_mode}; cd ${CURDIR};)
	-@$(foreach dir, ${DESTINY_alts}, rm -rf ${FINAL_DESTINY_ADJ}${dir};)
	@echo "Creating subdirectories in ${FINAL_DESTINY_ADJ}..."; $(foreach dir, ${DESTINY_alts}, mkdir -p ${FINAL_DESTINY_ADJ}${dir};)
	-@$(foreach dir, ${DESTINY_alts}, cp ${PROCESS_DIR}${dir}/adjoint/* ${LOG_MSG} ${FINAL_DESTINY_ADJ}${dir}/;)
	@-mkdir -p ${AD_MAKE_DIR}/../user/none/ad_${TOOL}/${phys_base}
	@cp ${AD_MAKE_DIR}/calc_user_ad.f90 ${AD_MAKE_DIR}/../user/none/ad_${TOOL}/${phys_base}/.


clean :
	-rm -f core
	-rm -rf ${PROCESS_DIR}

cleaner :
	-rm -f core
	-rm -rf ${PROCESS_DIR}
	-$(foreach dir, ${DESTINY_alts_mkdir}, rm -rf ${FINAL_DESTINY}${dir};)

showall:
	@echo DESTINY_alts: ${DESTINY_alts}
	@echo DESTINY_alts_mkdir: ${DESTINY_alts_mkdir}
	@echo SHEMAT_DIR: ${SHEMAT_DIR}
	@echo FINAL_DESTINY: ${FINAL_DESTINY}
	@echo PROCESS_DIR: ${PROCESS_DIR}
	@$(foreach dir, ${DESTINY_alts_mkdir}, echo ${FINAL_DESTINY}${dir};)

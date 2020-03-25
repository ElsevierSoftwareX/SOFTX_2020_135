# MIT License
#
# Copyright (c) 2020 SHEMAT-Suite
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Standard AD USER and PROPS
set(USERS_PROP "basc")
set(PROPS_USER "wells3d")

# AD Forward targets
#===================
add_custom_target(ad_exports
   COMMAND
   export forward_lst="${AD_DIRECTORIES}"; export users_lst="${USER_DIRECTORIES}";   export props_lst="${PROPS_DIRECTORIES}"; export props_user="${USERS_PROP}";export users_prop="${PROPS_USER}"
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mkAD
   COMMENT
     "export the propper ad bash variables"
)

add_custom_target(tap_tlm_users_prepare
   DEPENDS ad_exports
   COMMAND
   make TOOL=tap phys_base='${phys_base}' ad_mode='full' users_lst='${USER_DIRECTORIES}' props_lst='${USERS_PROP}' forward_lst='${AD_DIRECTORIES}' prepare
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mkAD
   COMMENT
      "prepare user functions for automatic differentiation with tapenade"
)


add_custom_target(tap_tlm_users_tlm
   DEPENDS tap_tlm_users_prepare
   COMMAND
   make TOOL=tap phys_base='${phys_base}' ad_mode='full' users_lst='${USER_DIRECTORIES}' props_lst='${USERS_PROP}' forward_lst='${AD_DIRECTORIES}' tlm
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mkAD
   COMMENT
      "execute user automatic differentiation with tapenade and copy the results"
)

add_custom_target(tap_tlm_props_prepare
   DEPENDS tap_tlm_users_tlm
   COMMAND
   make TOOL=tap phys_base='${phys_base}' ad_mode='full' users_lst='${PROPS_USER}' props_lst='${PROPS_DIRECTORIES}' forward_lst='${AD_DIRECTORIES}' prepare
   #echo "TOOL=tap phys_base='${phys_base}' ad_mode='full' users_lst='${USER}' props_lst='${USERS_PROP}' forward_lst='${AD_DIRECTORIES}' prepare"
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mkAD
   COMMENT
      "prepare props functions for automatic differentiation with tapenade"
)


add_custom_target(tap_tlm_props_tlm
   DEPENDS tap_tlm_props_prepare
   COMMAND
   make TOOL=tap phys_base='${phys_base}' ad_mode='full' users_lst='${PROPS_USER}' props_lst='${PROPS_DIRECTORIES}' forward_lst='${AD_DIRECTORIES}' tlm
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mkAD
   COMMENT
      "execute props automatic differentiation with tapenade and copy results"
)
add_custom_target(tap_tlm 
   DEPENDS tap_tlm_props_tlm
   COMMAND
   export forward_lst="${AD_DIRECTORIES}" users_lst="${USER_DIRECTORIES}" props_lst="${PROPS_DIRECTORIES}" props_user="${PROPS_USER}" users_prop=${USERS_PROP} phys_base="${phys_base}" && 
   ./mkAD/mv_tap_tlm && 
   rm user/none/g_tap/${phys_base}/calc_user_ftl.f90 &&
   rm -fr g_tap/${phys_base}/*_nodiff.f* g_tap/${phys_base}/daxpy_ftl.f* g_tap/${phys_base}/dcopy_ftl.f* g_tap/${phys_base}/dscal_ftl.f* props/*/g_tap/${phys_base}/*_nodiff*.* user/*/g_tap/${phys_base}/*_nodiff*.f90
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
   COMMENT "AD Differentiation with Tapenade"
   )

add_custom_target(tap_tlm_clean
   COMMAND
    rm -fr ${CMAKE_CURRENT_SOURCE_DIR}/props/*/g_tap/${phys_base} ${CMAKE_CURRENT_SOURCE_DIR}/user/*/g_tap/${phys_base} ${CMAKE_CURRENT_SOURCE_DIR}/g_tap/${phys_base}
    COMMENT "Cleaning Tapenade TLM directories"
    )


# AD Reverse targets
#===================
add_custom_target(tap_adm_users_prepare
   DEPENDS ad_exports
   COMMAND
   make TOOL=tap phys_base='${phys_base}' ad_mode='full' users_lst='${USER_DIRECTORIES}' props_lst='${USERS_PROP}' forward_lst='${AD_DIRECTORIES}' prepare_reverse
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mkAD
   COMMENT
      "prepare user functions for automatic differentiation with tapenade"
)


add_custom_target(tap_adm_users_adm
   DEPENDS tap_adm_users_prepare
   COMMAND
   make TOOL=tap phys_base='${phys_base}' ad_mode='full' users_lst='${USER_DIRECTORIES}' props_lst='${USERS_PROP}' forward_lst='${AD_DIRECTORIES}' adm
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mkAD
   COMMENT
      "execute user automatic differentiation with tapenade and copy the results"
)

add_custom_target(tap_adm_props_prepare
   DEPENDS tap_adm_users_adm
   COMMAND
   make TOOL=tap phys_base='${phys_base}' ad_mode='full' users_lst='${PROPS_USER}' props_lst='${PROPS_DIRECTORIES}' forward_lst='${AD_DIRECTORIES}' prepare_reverse
   #echo "TOOL=tap phys_base='${phys_base}' ad_mode='full' users_lst='${USER}' props_lst='${USERS_PROP}' forward_lst='${AD_DIRECTORIES}' prepare"
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mkAD
   COMMENT
      "prepare props functions for automatic differentiation with tapenade"
)


add_custom_target(tap_adm_props_adm
   DEPENDS tap_adm_props_prepare
   COMMAND
   make TOOL=tap phys_base='${phys_base}' ad_mode='full' users_lst='${PROPS_USER}' props_lst='${PROPS_DIRECTORIES}' forward_lst='${AD_DIRECTORIES}' adm
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mkAD
   COMMENT
      "execute props automatic differentiation with tapenade and copy results"
)
add_custom_target(tap_adm
   DEPENDS tap_adm_props_adm
   COMMAND
   export forward_lst="${AD_DIRECTORIES}" users_lst="${USER_DIRECTORIES}" props_lst="${PROPS_DIRECTORIES}" props_user="${PROPS_USER}" users_prop=${USERS_PROP} phys_base="${phys_base}" && 
   ./mkAD/mv_tap_adm #&& 
   #rm user/none/ad_tap/ad_calc_user.f90 &&
   #rm -fr ad_tap/*_nodiff.f* ad_tap/daxpy_ftl.f* ad_tap/dcopy_ftl.f* ad_tap/dscal_ftl.f* props/*/ad_tap/*_nodiff.* user/*/ad_tap/*_nodiff.f90
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
   COMMENT "AD Differentiation with Tapenade"
   )

add_custom_target(tap_adm_clean
   COMMAND
    rm -fr ${CMAKE_CURRENT_SOURCE_DIR}/props/*/ad_tap/${phys_base}\${CMAKE_CURRENT_SOURCE_DIR}/user/*/ad_tap/${phys_base} ${CMAKE_CURRENT_SOURCE_DIR}/ad_tap/${phys_base}
    COMMENT "Cleaning Tapenade TLM directories"
    )

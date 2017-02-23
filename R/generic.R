###################################################################
# planor R package
# Copyright INRA 2017
# INRA, UR1404, Research Unit MaIAGE
# F78350 Jouy-en-Josas, France.
#
# URL: http://www3.jouy.inra.fr/miaj/public/logiciels/planor/
#
# This file is part of planor R package.
# planor is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################

##---------------------------------------------------------------------------
# GENERIC FUNCTIONS
##---------------------------------------------------------------------------

# ---------------------------------------------------
# Generic "bind"
# ---------------------------------------------------
#  Generic function to bind two objects.
# ARGUMENTS
#  - x: first object
#  - y: second object
#  - ...: other arguments
# RETURN
#  An object where the input are binded.
# --------------------------------------
setGeneric("bind",
           function(x,y,...){ value <- standardGeneric("bind") },
           useAsDefault=FALSE)
## --------------------------------------------
## Generic "pick"
## --------------------------------------------
#  Generic function to extract one object from a complex object
# ARGUMENTS
#  - keys : an object from which to pick
#  - selection:  index vector to indicate the selection
# RETURN
# The selected object
## --------------------------------------------------------
setGeneric("pick",
           function(keys,...){
             value <- standardGeneric("pick")
         })
## --------------------------------------------
# Generic "planor.design" 
# --------------------------------------
#  Generic function to build a design
# ARGUMENTS
#   - key : an object from which the design will be built
# RETURN
# An object which  contains the design build from the input.
# --------------------------------------
setGeneric("planor.design",
           function(key, ...){
               value <- standardGeneric("planor.design")
           })
##---------------------------------------------------------------------------
setGeneric("getDesign",
           function(object, ...){
               value <- standardGeneric("getDesign")
           })

setGeneric("alias",
           function(object, model, ...){
             value <- standardGeneric("alias")
             })

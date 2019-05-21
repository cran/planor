###################################################################
# planor R package
# Copyright INRA 2019
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
.onAttach <- function(libname, pkgname)
  {
    packageStartupMessage("Loaded planor ", as.character(utils::packageVersion("planor")),"\n")
  }

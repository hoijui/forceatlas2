/*
SPDX-License-Identifier: GPL-3.0-only OR CDDL-1.0
Copyright 2008-2011 Gephi
Authors : Mathieu Jacomy <mathieu.jacomy@gmail.com>
Website : http://www.gephi.org

This file is part of Gephi.

DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS HEADER.

Copyright 2011 Gephi Consortium. All rights reserved.

The contents of this file are subject to the terms of either the GNU
General Public License Version 3 only ("GPL") or the Common
Development and Distribution License("CDDL") (collectively, the
"License"). You may not use this file except in compliance with the
License. You can obtain a copy of the License at
http://gephi.org/about/legal/license-notice/
or /cddl-1.0.txt and /gpl-3.0.txt. See the License for the
specific language governing permissions and limitations under the
License.  When distributing the software, include this License Header
Notice in each file and include the License files at
/cddl-1.0.txt and /gpl-3.0.txt. If applicable, add the following below the
License Header, with the fields enclosed by brackets [] replaced by
your own identifying information:
"Portions Copyrighted [year] [name of copyright owner]"

If you wish your version of this file to be governed by only the CDDL
or only the GPL Version 3, indicate your decision by adding
"[Contributor] elects to include this software in this distribution
under the [CDDL or GPL Version 3] license." If you do not indicate a
single choice of license, a recipient has the option to distribute
your version of this file under either the CDDL, the GPL Version 3 or
to extend the choice of license to its licensees as provided above.
However, if you add GPL Version 3 code and therefore, elected the GPL
Version 3 license, then the option applies only if the new code is
made subject to such option by the copyright holder.

Contributor(s):

Portions Copyrighted 2011 Gephi Consortium.
 */

package kco.forceatlas2;

/**
 * Data stored in Nodes and used by ForceAtlas2
 *
 * @author Mathieu Jacomy
 */
public class ForceAtlas2LayoutData2d implements ForceAtlas2LayoutData {

    private double dx;
    private double dy;
    private double old_dx;
    private double old_dy;
    private double mass = 1;

    @Override
    public double getDx() {
        return dx;
    }

    @Override
    public void setDx(final double dx) {
        this.dx = dx;
    }

    @Override
    public double getDy() {
        return dy;
    }

    @Override
    public void setDy(final double dy) {
        this.dy = dy;
    }

    @Override
    public double getDz() {
        return 0;
    }

    @Override
    public void setDz(final double dz) {
    }

    @Override
    public double getOld_dx() {
        return old_dx;
    }

    @Override
    public void setOld_dx(final double old_dx) {
        this.old_dx = old_dx;
    }

    @Override
    public double getOld_dy() {
        return old_dy;
    }

    @Override
    public void setOld_dy(final double old_dy) {
        this.old_dy = old_dy;
    }

    @Override
    public double getOld_dz() {
        return 0;
    }

    @Override
    public void setOld_dz(final double old_dz) {
    }

    @Override
    public double getMass() {
        return mass;
    }

    @Override
    public void setMass(final double mass) {
        this.mass = mass;
    }

    @Override
    public synchronized void augmentDx(double ddx) {
        this.dx += ddx;
    }

    @Override
    public synchronized void augmentDy(final double ddy) {
        this.dy += ddy;
    }

    @Override
    public void augmentDz(final double ddz) {
    }
}

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

import org.gephi.graph.api.Node;

/**
 * Generates the forces on demand, here are all the formulas for attraction and
 * repulsion.
 *
 * @author Mathieu Jacomy
 */
public class ForceFactory {

    public static final ForceFactory builder = new ForceFactory();

    private ForceFactory() {
    }

    public RepulsionForce buildRepulsion(final boolean adjustBySize, final double coefficient) {
        if (adjustBySize) {
            return new linRepulsion_antiCollision(coefficient);
        } else {
            return new linRepulsion(coefficient);
        }
    }

    public RepulsionForce getStrongGravity(final double coefficient) {
        return new strongGravity(coefficient);
    }

    public AttractionForce buildAttraction(
            final boolean logAttraction,
            final boolean distributedAttraction,
            final boolean adjustBySize,
            final double coefficient)
    {
        if (adjustBySize) {
            if (logAttraction) {
                if (distributedAttraction) {
                    return new logAttraction_degreeDistributed_antiCollision(coefficient);
                } else {
                    return new logAttraction_antiCollision(coefficient);
                }
            } else {
                if (distributedAttraction) {
                    return new linAttraction_degreeDistributed_antiCollision(coefficient);
                } else {
                    return new linAttraction_antiCollision(coefficient);
                }
            }
        } else {
            if (logAttraction) {
                if (distributedAttraction) {
                    return new logAttraction_degreeDistributed(coefficient);
                } else {
                    return new logAttraction(coefficient);
                }
            } else {
                if (distributedAttraction) {
                    return new linAttraction_massDistributed(coefficient);
                } else {
                    return new linAttraction(coefficient);
                }
            }
        }
    }

    public static abstract class AttractionForce {

        /**
         * Model for node-node attraction.
         *
         * @param n1
         * @param n2
         * @param e edge weight if needed
         */
        public abstract void apply(Node n1, Node n2, double e);
    }

    public static abstract class RepulsionForce {

        /**
         * Model for node-node repulsion
         * @param n1
         * @param n2
         */
        public abstract void apply(Node n1, Node n2);

        /**
         * Model for Barnes Hut approximation.
         *
         * @param n
         * @param r
         */
        public abstract void apply(Node n, Region r);

        /**
         * Model for gravitation (anti-repulsion).
         *
         * @param n the node to apply the gravitation force on
         * @param g the gravitation force/strength to apply ot the node
         */
        public abstract void apply(Node n, double g);

        /**
         * Model for node-node repulsion in quadtree (BH), do not update both n and o.
         *
         * @param n
         * @param o
         */
        public abstract void apply_BH(Node n, Node o);
    }

    /*
     * Repulsion force: Linear
     */
    private static class linRepulsion extends RepulsionForce {

        private final double coefficient;

        public linRepulsion(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2) {
            final ForceAtlas2LayoutData n1Layout = n1.getLayoutData();
            final ForceAtlas2LayoutData n2Layout = n2.getLayoutData();

            // Get the distance

            final double xDist = n1.x() - n2.x();
            final double yDist = n1.y() - n2.y();
            final double zDist = n1.z() - n2.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = coefficient * n1Layout.getMass() * n2Layout.getMass() / distance / distance;

                n1Layout.setDx(n1Layout.getDx() + xDist * factor);
                n1Layout.setDy(n1Layout.getDy() + yDist * factor);
                n1Layout.setDz(n1Layout.getDz() + zDist * factor);

                n2Layout.setDx(n2Layout.getDx() - xDist * factor);
                n2Layout.setDy(n2Layout.getDy() - yDist * factor);
                n2Layout.setDz(n2Layout.getDz() - zDist * factor);
            }
        }

        @Override
        public void apply(final Node n, final Region r) {
            final ForceAtlas2LayoutData nLayout = n.getLayoutData();

            // Get the distance
            final double xDist = n.x() - r.getMassCenterX();
            final double yDist = n.y() - r.getMassCenterY();
            final double zDist = n.z() - r.getMassCenterZ();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = coefficient * nLayout.getMass() * r.getMass() / distance / distance;

                nLayout.setDx(nLayout.getDx() + xDist * factor);
                nLayout.setDy(nLayout.getDy() + yDist * factor);
                nLayout.setDz(nLayout.getDz() + zDist * factor);
            }
        }

        @Override
        public void apply(final Node n, final double g) {
            final ForceAtlas2LayoutData nLayout = n.getLayoutData();

            // Get the distance
            final double xDist = n.x();
            final double yDist = n.y();
            final double zDist = n.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = coefficient * nLayout.getMass() * g / distance;

                nLayout.setDx(nLayout.getDx() - xDist * factor);
                nLayout.setDy(nLayout.getDy() - yDist * factor);
                nLayout.setDz(nLayout.getDz() - zDist * factor);
            }
        }

        @Override
        public void apply_BH(final Node n, final Node o) {
            final ForceAtlas2LayoutData nLayout = n.getLayoutData();
            final ForceAtlas2LayoutData oLayout = o.getLayoutData();

            // Get the distance
            final double xDist = n.x() - o.x();
            final double yDist = n.y() - o.y();
            final double zDist = n.z() - o.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = coefficient * nLayout.getMass() * oLayout.getMass() / distance / distance;

                nLayout.setDx(nLayout.getDx() + xDist * factor);
                nLayout.setDy(nLayout.getDy() + yDist * factor);
                nLayout.setDz(nLayout.getDz() + zDist * factor);
            }
        }
    }

    /*
     * Repulsion force: Strong Gravity (as a Repulsion Force because it is easier)
     */
    private static class linRepulsion_antiCollision extends RepulsionForce {

        private final double coefficient;

        public linRepulsion_antiCollision(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2) {
            final ForceAtlas2LayoutData n1Layout = n1.getLayoutData();
            final ForceAtlas2LayoutData n2Layout = n2.getLayoutData();

            // Get the distance
            final double xDist = n1.x() - n2.x();
            final double yDist = n1.y() - n2.y();
            final double zDist = n1.z() - n2.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist) - n1.size() - n2.size();

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = coefficient * n1Layout.getMass() * n2Layout.getMass() / distance / distance;

                n1Layout.setDx(n1Layout.getDx() + xDist * factor);
                n1Layout.setDy(n1Layout.getDy() + yDist * factor);
                n1Layout.setDz(n1Layout.getDz() + zDist * factor);

                n2Layout.setDx(n2Layout.getDx() - xDist * factor);
                n2Layout.setDy(n2Layout.getDy() - yDist * factor);
                n2Layout.setDz(n2Layout.getDz() - zDist * factor);

            } else if (distance < 0) {
                final double factor = 100 * coefficient * n1Layout.getMass() * n2Layout.getMass();

                n1Layout.setDx(n1Layout.getDx() + xDist * factor);
                n1Layout.setDy(n1Layout.getDy() + yDist * factor);
                n1Layout.setDz(n1Layout.getDz() + zDist * factor);

                n2Layout.setDx(n2Layout.getDx() - xDist * factor);
                n2Layout.setDy(n2Layout.getDy() - yDist * factor);
                n2Layout.setDz(n2Layout.getDz() - zDist * factor);
            }
        }

        @Override
        public void apply(final Node n, final Region r) {
            final ForceAtlas2LayoutData nLayout = n.getLayoutData();

            // Get the distance
            final double xDist = n.x() - r.getMassCenterX();
            final double yDist = n.y() - r.getMassCenterY();
            final double zDist = n.z() - r.getMassCenterZ();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = coefficient * nLayout.getMass() * r.getMass() / distance / distance;

                nLayout.setDx(nLayout.getDx() + xDist * factor);
                nLayout.setDy(nLayout.getDy() + yDist * factor);
                nLayout.setDz(nLayout.getDz() + zDist * factor);
            } else if (distance < 0) {
                final double factor = -coefficient * nLayout.getMass() * r.getMass() / distance;

                nLayout.setDx(nLayout.getDx() + xDist * factor);
                nLayout.setDy(nLayout.getDy() + yDist * factor);
                nLayout.setDz(nLayout.getDz() + zDist * factor);
            }
        }

        @Override
        public void apply(final Node n, final double g) {
            final ForceAtlas2LayoutData nLayout = n.getLayoutData();

            // Get the distance
            final double xDist = n.x();
            final double yDist = n.y();
            final double zDist = n.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = coefficient * nLayout.getMass() * g / distance;

                nLayout.setDx(nLayout.getDx() - xDist * factor);
                nLayout.setDy(nLayout.getDy() - yDist * factor);
                nLayout.setDz(nLayout.getDz() - zDist * factor);
            }
        }

        @Override
        public void apply_BH(final Node n, final Node o) {
            final ForceAtlas2LayoutData nLayout = n.getLayoutData();
            final ForceAtlas2LayoutData oLayout = o.getLayoutData();
            // Get the distance
            final double xDist = n.x() - o.x();
            final double yDist = n.y() - o.y();
            final double zDist = n.z() - o.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist) - n.size() - o.size();

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = coefficient * nLayout.getMass() * oLayout.getMass() / distance / distance;

                nLayout.setDx(nLayout.getDx() + xDist * factor);
                nLayout.setDy(nLayout.getDy() + yDist * factor);
                nLayout.setDz(nLayout.getDz() + zDist * factor);
            }
        }
    }

    private static class strongGravity extends RepulsionForce {

        private final double coefficient;

        public strongGravity(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2) {
            // Not Relevant
        }

        @Override
        public void apply(final Node n, final Region r) {
            // Not Relevant
        }

        @Override
        public void apply(final Node n, final double g) {
            final ForceAtlas2LayoutData nLayout = n.getLayoutData();

            // Get the distance
            final double xDist = n.x();
            final double yDist = n.y();
            final double zDist = n.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = coefficient * nLayout.getMass() * g;

                nLayout.setDx(nLayout.getDx() - xDist * factor);
                nLayout.setDy(nLayout.getDy() - yDist * factor);
                nLayout.setDz(nLayout.getDz() - zDist * factor);
            }
        }

        @Override
        public void apply_BH(final Node n, final Node o) {
            // Not Relevant
        }
    }

    /*
     * Attraction force: Linear
     */
    private static class linAttraction extends AttractionForce {

        private final double coefficient;

        public linAttraction(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2, final double e) {
            final ForceAtlas2LayoutData n1Layout = n1.getLayoutData();
            final ForceAtlas2LayoutData n2Layout = n2.getLayoutData();

            // Get the distance
            final double xDist = n1.x() - n2.x();
            final double yDist = n1.y() - n2.y();
            final double zDist = n1.z() - n2.z();

            // NB: factor = force / distance
            final double factor = -coefficient * e;

            // if (n1.getId().equals("10812")) {
            //     System.out.println("before (" + n2.getId() + ") " + n1Layout.getDx() + " " + n1Layout.getDy());
            // }
            n1Layout.augmentDx(xDist * factor);
            n1Layout.augmentDy(yDist * factor);
            n1Layout.augmentDz(zDist * factor);
            // if (n1.getId().equals("10812")) {
            //     System.out.println("after(" + n2.getId() +  ") " + n1Layout.getDx() + " " + n1Layout.getDy());
            // }

            // if (n2.getId().equals("10812")) {
            //     System.out.println("before (" + n1.getId() + ") " + n2Layout.getDx() + " " + n2Layout.getDy());
            // }
            n2Layout.augmentDx(- xDist * factor);
            n2Layout.augmentDy(- yDist * factor);
            n2Layout.augmentDz(- zDist * factor);
            // if (n2.getId().equals("10812")) {
            //     System.out.println("after (" + n1.getId() + ") " + n2Layout.getDx() + " " + n2Layout.getDy());
            // }
        }
    }

    /*
     * Attraction force: Linear, distributed by mass (typically, degree)
     */
    private static class linAttraction_massDistributed extends AttractionForce {

        private final double coefficient;

        public linAttraction_massDistributed(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2, final double e) {
            final ForceAtlas2LayoutData n1Layout = n1.getLayoutData();
            final ForceAtlas2LayoutData n2Layout = n2.getLayoutData();

            // Get the distance
            final double xDist = n1.x() - n2.x();
            final double yDist = n1.y() - n2.y();
            final double zDist = n1.z() - n2.z();

            // NB: factor = force / distance
            final double factor = -coefficient * e / n1Layout.getMass();

            n1Layout.augmentDx(xDist * factor);
            n1Layout.augmentDy(yDist * factor);
            n1Layout.augmentDz(zDist * factor);

            n2Layout.augmentDx(- xDist * factor);
            n2Layout.augmentDy(- yDist * factor);
            n2Layout.augmentDz(- zDist * factor);
        }
    }

    /*
     * Attraction force: Logarithmic
     */
    private static class logAttraction extends AttractionForce {

        private final double coefficient;

        public logAttraction(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2, final double e) {
            final ForceAtlas2LayoutData n1Layout = n1.getLayoutData();
            final ForceAtlas2LayoutData n2Layout = n2.getLayoutData();

            // Get the distance
            final double xDist = n1.x() - n2.x();
            final double yDist = n1.y() - n2.y();
            final double zDist = n1.z() - n2.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = -coefficient * e * Math.log(1 + distance) / distance;

                n1Layout.augmentDx(xDist * factor);
                n1Layout.augmentDy(yDist * factor);
                n1Layout.augmentDz(zDist * factor);

                n2Layout.augmentDx(- xDist * factor);
                n2Layout.augmentDy(- yDist * factor);
                n2Layout.augmentDz(- zDist * factor);
            }
        }
    }

    /*
     * Attraction force: Linear, distributed by Degree
     */
    private static class logAttraction_degreeDistributed extends AttractionForce {

        private final double coefficient;

        public logAttraction_degreeDistributed(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2, final double e) {
            final ForceAtlas2LayoutData n1Layout = n1.getLayoutData();
            final ForceAtlas2LayoutData n2Layout = n2.getLayoutData();

            // Get the distance
            final double xDist = n1.x() - n2.x();
            final double yDist = n1.y() - n2.y();
            final double zDist = n1.z() - n2.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = -coefficient * e * Math.log(1 + distance) / distance / n1Layout.getMass();

                n1Layout.augmentDx(xDist * factor);
                n1Layout.augmentDy(yDist * factor);
                n1Layout.augmentDz(zDist * factor);

                n2Layout.augmentDx(- xDist * factor);
                n2Layout.augmentDy(- yDist * factor);
                n2Layout.augmentDz(- zDist * factor);
            }
        }
    }

    /*
     * Attraction force: Linear, with Anti-Collision
     */
    private static class linAttraction_antiCollision extends AttractionForce {

        private final double coefficient;

        public linAttraction_antiCollision(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2, final double e) {
            final ForceAtlas2LayoutData n1Layout = n1.getLayoutData();
            final ForceAtlas2LayoutData n2Layout = n2.getLayoutData();

            // Get the distance
            final double xDist = n1.x() - n2.x();
            final double yDist = n1.y() - n2.y();
            final double zDist = n1.z() - n2.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist) - n1.size() - n2.size();

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = -coefficient * e;

                n1Layout.augmentDx(xDist * factor);
                n1Layout.augmentDy(yDist * factor);
                n1Layout.augmentDz(zDist * factor);

                n2Layout.augmentDx(- xDist * factor);
                n2Layout.augmentDy(- yDist * factor);
                n2Layout.augmentDz(- zDist * factor);
            }
        }
    }

    /*
     * Attraction force: Linear, distributed by Degree, with Anti-Collision
     */
    private static class linAttraction_degreeDistributed_antiCollision extends AttractionForce {

        private final double coefficient;

        public linAttraction_degreeDistributed_antiCollision(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2, final double e) {
            final ForceAtlas2LayoutData n1Layout = n1.getLayoutData();
            final ForceAtlas2LayoutData n2Layout = n2.getLayoutData();

            // Get the distance
            final double xDist = n1.x() - n2.x();
            final double yDist = n1.y() - n2.y();
            final double zDist = n1.z() - n2.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist) - n1.size() - n2.size();

            if (distance > 0) {
                // NB: factor = force / distance
                final double factor = -coefficient * e / n1Layout.getMass();

                n1Layout.augmentDx(xDist * factor);
                n1Layout.augmentDy(yDist * factor);
                n1Layout.augmentDz(zDist * factor);

                n2Layout.augmentDx(- xDist * factor);
                n2Layout.augmentDy(- yDist * factor);
                n2Layout.augmentDz(- zDist * factor);
            }
        }
    }

    /*
     * Attraction force: Logarithmic, with Anti-Collision
     */
    private static class logAttraction_antiCollision extends AttractionForce {

        private final double coefficient;

        public logAttraction_antiCollision(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2, final double e) {
            final ForceAtlas2LayoutData n1Layout = n1.getLayoutData();
            final ForceAtlas2LayoutData n2Layout = n2.getLayoutData();

            // Get the distance
            final double xDist = n1.x() - n2.x();
            final double yDist = n1.y() - n2.y();
            final double zDist = n1.z() - n2.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist) - n1.size() - n2.size();

            if (distance > 0) {

                // NB: factor = force / distance
                final double factor = -coefficient * e * Math.log(1 + distance) / distance;

                n1Layout.augmentDx(xDist * factor);
                n1Layout.augmentDy(yDist * factor);
                n1Layout.augmentDz(zDist * factor);

                n2Layout.augmentDx(- xDist * factor);
                n2Layout.augmentDy(- yDist * factor);
                n2Layout.augmentDz(- zDist * factor);
            }
        }
    }

    /*
     * Attraction force: Linear, distributed by Degree, with Anti-Collision
     */
    private static class logAttraction_degreeDistributed_antiCollision extends AttractionForce {

        private final double coefficient;

        public logAttraction_degreeDistributed_antiCollision(final double c) {
            coefficient = c;
        }

        @Override
        public void apply(final Node n1, final Node n2, final double e) {
            final ForceAtlas2LayoutData n1Layout = n1.getLayoutData();
            final ForceAtlas2LayoutData n2Layout = n2.getLayoutData();

            // Get the distance
            final double xDist = n1.x() - n2.x();
            final double yDist = n1.y() - n2.y();
            final double zDist = n1.z() - n2.z();
            final double distance = Math.sqrt(xDist * xDist + yDist * yDist + zDist * zDist) - n1.size() - n2.size();

            if (distance > 0) {

                // NB: factor = force / distance
                final double factor = -coefficient * e * Math.log(1 + distance) / distance / n1Layout.getMass();

                n1Layout.augmentDx(xDist * factor);
                n1Layout.augmentDy(yDist * factor);
                n1Layout.augmentDz(zDist * factor);

                n2Layout.augmentDx(- xDist * factor);
                n2Layout.augmentDy(- yDist * factor);
                n2Layout.augmentDz(- zDist * factor);
            }
        }
    }
}

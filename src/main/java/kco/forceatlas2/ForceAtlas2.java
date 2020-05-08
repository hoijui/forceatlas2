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

import org.gephi.graph.api.Edge;
import org.gephi.graph.api.Graph;
import org.gephi.graph.api.GraphModel;
import org.gephi.graph.api.Interval;
import org.gephi.graph.api.Node;
import org.gephi.layout.spi.Layout;
import org.gephi.layout.spi.LayoutBuilder;
import org.gephi.layout.spi.LayoutProperty;
import org.openide.util.Exceptions;
import org.openide.util.NbBundle;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * ForceAtlas 2 Layout, manages each step of the computations.
 *
 * @author Mathieu Jacomy
 * @author Joshua Gould
 */
public class ForceAtlas2 implements Layout {

    private final ForceAtlas2Builder layoutBuilder;
    private double outboundAttCompensation = 1;
    private GraphModel graphModel;
    private Graph graph;
    private double edgeWeightInfluence;
    private double jitterTolerance;
    private double scalingRatio;
    private double gravity;
    private double speed;
    private double speedEfficiency;
    private boolean outboundAttractionDistribution;
    private boolean adjustSizes;
    private boolean barnesHutOptimize;
    private double barnesHutTheta;
    private boolean linLogMode;
    private boolean strongGravityMode;
    private int threadCount;
    private int currentThreadCount;
    private int stepCount;
    private int updateBarnesHutIter = 1;
    private Region rootRegion;
    private ExecutorService pool;
    private boolean updateCenter = true;
    private Node[] nodes;
    private Edge[] edges;
    private int barnesHutSplits = -1;
    private double distance;
    private final boolean is3d;
    private final boolean useAltSpeed;

    public ForceAtlas2(final ForceAtlas2Builder layoutBuilder, final boolean is3d, final boolean useAltSpeed) {
        this.layoutBuilder = layoutBuilder;
        this.is3d = is3d;
        this.useAltSpeed = useAltSpeed;
        this.threadCount = Runtime.getRuntime().availableProcessors();
    }

    private static double getEdgeWeight(final Edge edge, final boolean isDynamicWeight, final Interval interval) {
        if (isDynamicWeight) {
            return edge.getWeight(interval);
        } else {
            return edge.getWeight();
        }
    }

    private static void waitForFutures(final List<Future<?>> futures) {
        for (Future<?> f : futures) {
            try {
                f.get();
            } catch (Exception x) {
                x.printStackTrace();
            }
        }
    }

    @Override
    public void initAlgo() {
//        AbstractLayout.ensureSafeLayoutNodePositions(graphModel);
        stepCount = 0;
        speed = 1.;
        speedEfficiency = 1.;
        graph = graphModel.getGraphVisible();
        nodes = graph.getNodes().toArray();
        edges = graph.getEdges().toArray();
        pool = Executors.newFixedThreadPool(threadCount);
        if (this.barnesHutSplits == -1) {
            this.barnesHutSplits = (int) Math.floor(Math.log(this.threadCount) / Math.log(is3d ? 8.0 : 4.0) + 0.02) + 1;
        }
        // Initialise layout data
        for (Node n : nodes) {
            if (n.getLayoutData() == null) {
                ForceAtlas2LayoutData nLayout = is3d ? new ForceAtlas2LayoutData3d() : new ForceAtlas2LayoutData2d();
                n.setLayoutData(nLayout);
            }
            ForceAtlas2LayoutData nLayout = n.getLayoutData();
            nLayout.setMass(1 + graph.getDegree(n));
            nLayout.setOld_dx(0);
            nLayout.setOld_dy(0);
            nLayout.setOld_dz(0);
            nLayout.setDx(0);
            nLayout.setDy(0);
            nLayout.setDz(0);
        }

        currentThreadCount = threadCount;
    }

    private void barnesHutRegions() {

        if (isBarnesHutOptimize()) {
            if (stepCount % updateBarnesHutIter == 0) {
                rootRegion = new Region(nodes, is3d);
                List<Region> regions = new ArrayList<>();
                regions.add(rootRegion);
                for (int splitIndex = 0; splitIndex < barnesHutSplits; splitIndex++) {
                    List<Future<?>> futures = new ArrayList<>();
                    for (Region r : regions) {
                        futures.add(pool.submit(new BarnesHutBuildSubRegionTask(Collections.singletonList(r), false)));
                    }
                    waitForFutures(futures);
                    List<Region> newRegions = new ArrayList<>();
                    for (Region r : regions) {
                        newRegions.addAll(r.getSubRegions());
                    }
                    regions = newRegions;
                }

                List<Future<?>> futures = new ArrayList<>();
                for (int t = currentThreadCount; t > 0; t--) {
                    int from = (int) Math.floor((double) regions.size() * (t - 1) / currentThreadCount);
                    int to = (int) Math.floor((double) regions.size() * t / currentThreadCount);
                    futures.add(pool.submit(new BarnesHutBuildSubRegionTask(regions.subList(from, to), true)));
                }
                waitForFutures(futures);

            } else if (updateCenter) {
                rootRegion = new Region(nodes, is3d);
                List<Region> regions = new ArrayList<>();
                regions.add(rootRegion);
                for (int splitIndex = 0; splitIndex < barnesHutSplits; splitIndex++) {
                    List<Future<?>> futures = new ArrayList<>();
                    for (Region r : regions) {
                        futures.add(pool.submit(new BarnesHutUpdateCenterTask(Collections.singletonList(r), false)));
                    }
                    waitForFutures(futures);
                    List<Region> newRegions = new ArrayList<>();
                    for (Region r : regions) {
                        newRegions.addAll(r.getSubRegions());
                    }
                    regions = newRegions;
                }
                List<Future<?>> futures = new ArrayList<>();
                for (int t = currentThreadCount; t > 0; t--) {
                    int from = (int) Math.floor((double) regions.size() * (t - 1) / currentThreadCount);
                    int to = (int) Math.floor((double) regions.size() * t / currentThreadCount);
                    futures.add(pool.submit(new BarnesHutUpdateCenterTask(regions.subList(from, to), true)));
                }
                waitForFutures(futures);
            }
        }
    }

    private void repulsionAndGravity() {

        // Repulsion (and gravity)
        // NB: Multi-threaded
        ForceFactory.RepulsionForce Repulsion = ForceFactory.builder.buildRepulsion(isAdjustSizes(), getScalingRatio());

        List<Future<?>> futures = new ArrayList<>();
        for (int t = currentThreadCount; t > 0; t--) {
            int from = (int) Math.floor((double) nodes.length * (t - 1) / currentThreadCount);
            int to = (int) Math.floor((double) nodes.length * t / currentThreadCount);
            futures.add(pool.submit(new NodesThread(nodes, from, to, isBarnesHutOptimize(), getBarnesHutTheta(), getGravity(), (isStrongGravityMode()) ? (ForceFactory.builder.getStrongGravity(getScalingRatio())) : (Repulsion), getScalingRatio(), rootRegion, Repulsion)));
        }
        waitForFutures(futures);
    }

    private void outboundAttractionDistribution() {

        // If outboundAttractionDistribution active, compensate.
        if (isOutboundAttractionDistribution()) {
            outboundAttCompensation = 0;
            for (Node n : nodes) {
                ForceAtlas2LayoutData nLayout = n.getLayoutData();
                outboundAttCompensation += nLayout.getMass();
            }
            outboundAttCompensation /= nodes.length;
        }
    }

    private double applyForces() {

        List<Future<Double>> futures = new ArrayList<>();
        boolean adjustSizes = isAdjustSizes();
        List<Node> nodesList = Arrays.asList(nodes);
        for (int t = currentThreadCount; t > 0; t--) {
            int from = (int) Math.floor((double) nodes.length * (t - 1) / currentThreadCount);
            int to = (int) Math.floor((double) nodes.length * t / currentThreadCount);
            futures.add(pool.submit(new ApplyForcesTask(nodesList.subList(from, to), adjustSizes, speed, useAltSpeed)));
        }

        double distance = 0;
        try {
            for (Future<Double> f : futures) {
                distance += f.get();
            }
        } catch (Exception x) {
            x.printStackTrace();
        }

        return distance;
    }

    private void attraction(final boolean isDynamicWeight, final Interval interval) {
        List<Future<?>> futures = new ArrayList<>();

        final ForceFactory.AttractionForce Attraction = ForceFactory.builder.buildAttraction(isLinLogMode(), isOutboundAttractionDistribution(), isAdjustSizes(), 1 * ((isOutboundAttractionDistribution()) ? (outboundAttCompensation) : (1)));
        int taskCount = currentThreadCount;
        List<Edge> edgeList = Arrays.asList(edges);
        final Double edgeWeightInfluence = getEdgeWeightInfluence();
        for (int t = taskCount; t > 0; t--) {
            int from = (int) Math.floor((double) edges.length * (t - 1) / taskCount);
            int to = (int) Math.floor((double) edges.length * t / taskCount);
            final List<Edge> subList = edgeList.subList(from, to);
            Runnable task = new AttractionTask(subList, Attraction, isDynamicWeight, interval, edgeWeightInfluence);
            futures.add(pool.submit(task));
        }
        waitForFutures(futures);
    }

    private void speed() {

        List<Future<Double[]>> futures = new ArrayList<>();
        List<Node> nodesList = Arrays.asList(nodes);
        for (int t = currentThreadCount; t > 0; t--) {
            int from = (int) Math.floor((double) nodes.length * (t - 1) / currentThreadCount);
            int to = (int) Math.floor((double) nodes.length * t / currentThreadCount);
            futures.add(pool.submit(new SpeedTask(nodesList.subList(from, to))));
        }
        double totalSwinging = 0d;  // How much irregular movement
        double totalEffectiveTraction = 0d;  // Hom much useful movement
        try {
            for (final Future<Double[]> f : futures) {
                Double[] result = f.get();
                totalSwinging += result[0];
                totalEffectiveTraction += result[1];
            }
        } catch (InterruptedException | ExecutionException exc) {
            exc.printStackTrace();
        }

        // We want that swingingMovement < tolerance * convergenceMovement

        // Optimize jitter tolerance
        // The 'right' jitter tolerance for this network. Bigger networks need more tolerance. Denser networks need less tolerance. Totally empiric.
        double estimatedOptimalJitterTolerance = 0.05 * Math.sqrt(nodes.length);
        double minJT = Math.sqrt(estimatedOptimalJitterTolerance);
        double maxJT = 10;
        double jt = jitterTolerance * Math.max(minJT, Math.min(maxJT, estimatedOptimalJitterTolerance * totalEffectiveTraction / Math.pow(nodes.length, 2)));

        double minSpeedEfficiency = 0.05;

        // Protection against erratic behavior
        if (totalSwinging / totalEffectiveTraction > 2.0) {
            if (speedEfficiency > minSpeedEfficiency) {
                speedEfficiency *= 0.5;
            }
            jt = Math.max(jt, jitterTolerance);
        }

        double targetSpeed = jt * speedEfficiency * totalEffectiveTraction / totalSwinging;

        // Speed efficiency is how the speed really corresponds to the swinging vs. convergence tradeoff
        // We adjust it slowly and carefully
        if (totalSwinging > jt * totalEffectiveTraction) {
            if (speedEfficiency > minSpeedEfficiency) {
                speedEfficiency *= 0.7;
            }
        } else if (speed < 1000) {
            speedEfficiency *= 1.3;
        }

        // But the speed should not rise too much too quickly, since it would make the convergence drop dramatically.
        double maxRise = 0.5;   // Max rise: 50%
        speed = speed + Math.min(targetSpeed - speed, maxRise * speed);
    }

    private void initLayoutData() {

        List<Future<?>> futures = new ArrayList<>();
        List<Node> nodesList = Arrays.asList(nodes);
        for (int t = currentThreadCount; t > 0; t--) {
            int from = (int) Math.floor((double) nodes.length * (t - 1) / currentThreadCount);
            int to = (int) Math.floor((double) nodes.length * t / currentThreadCount);
            futures.add(pool.submit(new InitLayoutTask(nodesList.subList(from, to), graph)));
        }
        waitForFutures(futures);
    }

    @Override
    public void goAlgo() {

        boolean isDynamicWeight = graphModel.getEdgeTable().getColumn("weight").isDynamic();
        Interval interval = graph.getView().getTimeInterval();

        // Initialise layout data

        initLayoutData();

        // If Barnes Hut active, initialize root region
        barnesHutRegions();

        outboundAttractionDistribution();

        repulsionAndGravity();

        // Attraction
        attraction(isDynamicWeight, interval);

        // Auto adjust speed
        speed();

        // Apply forces
        distance = applyForces();
        stepCount++;
    }

    public double getDistance() {
        return distance;
    }

    @Override
    public boolean canAlgo() {
        return graphModel != null;
    }

    @Override
    public void endAlgo() {
//        for (Node n : graph.getNodes()) {
//            n.setLayoutData(null);
//        }

        pool.shutdown();
    }

    @Override
    public void resetPropertiesValues() {
        int nodesCount = 0;

        if (graphModel != null) {
            nodesCount = graphModel.getGraphVisible().getNodeCount();
        }

        // Tuning
        if (nodesCount >= 100) {
            setScalingRatio(2.0);
        } else {
            setScalingRatio(10.0);
        }
        setStrongGravityMode(false);
        setGravity(1.);

        // Behavior
        setOutboundAttractionDistribution(false);
        setLinLogMode(false);
        setAdjustSizes(false);
        setEdgeWeightInfluence(1.);

        // Performance
        setJitterTolerance(1d);
        setBarnesHutOptimize(nodesCount >= 1000);
        setBarnesHutTheta(1.2);
        setThreadsCount(Runtime.getRuntime().availableProcessors());
    }

    @Override
    public LayoutBuilder getBuilder() {
        return layoutBuilder;
    }

    @Override
    public void setGraphModel(final GraphModel graphModel) {
        this.graphModel = graphModel;
        // Trick: reset here to take the profile of the graph in account for default values
        resetPropertiesValues();
    }

    public Double getBarnesHutTheta() {
        return barnesHutTheta;
    }

    public void setBarnesHutTheta(final Double barnesHutTheta) {
        this.barnesHutTheta = barnesHutTheta;
    }

    public Double getEdgeWeightInfluence() {
        return edgeWeightInfluence;
    }

    public void setEdgeWeightInfluence(final Double edgeWeightInfluence) {
        this.edgeWeightInfluence = edgeWeightInfluence;
    }

    public Double getJitterTolerance() {
        return jitterTolerance;
    }

    public void setJitterTolerance(final Double jitterTolerance) {
        this.jitterTolerance = jitterTolerance;
    }

    public Boolean isLinLogMode() {
        return linLogMode;
    }

    public void setLinLogMode(final Boolean linLogMode) {
        this.linLogMode = linLogMode;
    }

    public Double getScalingRatio() {
        return scalingRatio;
    }

    public void setScalingRatio(final Double scalingRatio) {
        this.scalingRatio = scalingRatio;
    }

    public Boolean isStrongGravityMode() {
        return strongGravityMode;
    }

    public void setStrongGravityMode(final Boolean strongGravityMode) {
        this.strongGravityMode = strongGravityMode;
    }

    public Double getGravity() {
        return gravity;
    }

    public void setGravity(final Double gravity) {
        this.gravity = gravity;
    }

    public Integer getThreadsCount() {
        return threadCount;
    }

    public void setThreadsCount(final Integer threadCount) {
        this.threadCount = Math.max(1, threadCount);
    }

    public Boolean isOutboundAttractionDistribution() {
        return outboundAttractionDistribution;
    }

    public void setOutboundAttractionDistribution(final Boolean outboundAttractionDistribution) {
        this.outboundAttractionDistribution = outboundAttractionDistribution;
    }

    public Boolean isAdjustSizes() {
        return adjustSizes;
    }

    public void setAdjustSizes(final Boolean adjustSizes) {
        this.adjustSizes = adjustSizes;
    }

    public Boolean isBarnesHutOptimize() {
        return barnesHutOptimize;
    }

    public void setBarnesHutOptimize(Boolean barnesHutOptimize) {
        this.barnesHutOptimize = barnesHutOptimize;
    }

    public int getUpdateBarnesHutIter() {
        return updateBarnesHutIter;
    }

    public void setUpdateBarnesHutIter(final int updateBarnesHutIter) {
        this.updateBarnesHutIter = updateBarnesHutIter;
    }

    public void setUpdateCenter(final boolean updateCenter) {
        this.updateCenter = updateCenter;
    }

    public void setBarnesHutSplits(final int barnesHutSplits) {
        this.barnesHutSplits = barnesHutSplits;
    }

    @Override
    public LayoutProperty[] getProperties() {
        List<LayoutProperty> properties = new ArrayList<>();
        final String FORCE_ATLAS2_TUNING = NbBundle.getMessage(getClass(), "ForceAtlas2.tuning");
        final String FORCE_ATLAS2_BEHAVIOR = NbBundle.getMessage(getClass(), "ForceAtlas2.behavior");
        final String FORCE_ATLAS2_PERFORMANCE = NbBundle.getMessage(getClass(), "ForceAtlas2.performance");
        final String FORCE_ATLAS2_THREADS = NbBundle.getMessage(getClass(), "ForceAtlas2.threads");

        try {
            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.scalingRatio.name"),
                    FORCE_ATLAS2_TUNING,
                    "ForceAtlas2.scalingRatio.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.scalingRatio.desc"),
                    "getScalingRatio", "setScalingRatio"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.strongGravityMode.name"),
                    FORCE_ATLAS2_TUNING,
                    "ForceAtlas2.strongGravityMode.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.strongGravityMode.desc"),
                    "isStrongGravityMode", "setStrongGravityMode"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.gravity.name"),
                    FORCE_ATLAS2_TUNING,
                    "ForceAtlas2.gravity.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.gravity.desc"),
                    "getGravity", "setGravity"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.distributedAttraction.name"),
                    FORCE_ATLAS2_BEHAVIOR,
                    "ForceAtlas2.distributedAttraction.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.distributedAttraction.desc"),
                    "isOutboundAttractionDistribution", "setOutboundAttractionDistribution"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.linLogMode.name"),
                    FORCE_ATLAS2_BEHAVIOR,
                    "ForceAtlas2.linLogMode.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.linLogMode.desc"),
                    "isLinLogMode", "setLinLogMode"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.adjustSizes.name"),
                    FORCE_ATLAS2_BEHAVIOR,
                    "ForceAtlas2.adjustSizes.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.adjustSizes.desc"),
                    "isAdjustSizes", "setAdjustSizes"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.edgeWeightInfluence.name"),
                    FORCE_ATLAS2_BEHAVIOR,
                    "ForceAtlas2.edgeWeightInfluence.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.edgeWeightInfluence.desc"),
                    "getEdgeWeightInfluence", "setEdgeWeightInfluence"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.jitterTolerance.name"),
                    FORCE_ATLAS2_PERFORMANCE,
                    "ForceAtlas2.jitterTolerance.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.jitterTolerance.desc"),
                    "getJitterTolerance", "setJitterTolerance"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.barnesHutOptimization.name"),
                    FORCE_ATLAS2_PERFORMANCE,
                    "ForceAtlas2.barnesHutOptimization.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.barnesHutOptimization.desc"),
                    "isBarnesHutOptimize", "setBarnesHutOptimize"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.barnesHutTheta.name"),
                    FORCE_ATLAS2_PERFORMANCE,
                    "ForceAtlas2.barnesHutTheta.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.barnesHutTheta.desc"),
                    "getBarnesHutTheta", "setBarnesHutTheta"));

            properties.add(LayoutProperty.createProperty(
                    this, Integer.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas2.threads.name"),
                    FORCE_ATLAS2_THREADS,
                    "ForceAtlas2.threads.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas2.threads.desc"),
                    "getThreadsCount", "setThreadsCount"));
        } catch (NoSuchMethodException exc) {
            Exceptions.printStackTrace(exc);
        }

        return properties.toArray(new LayoutProperty[0]);
    }

    @Override
    public String toString() {
        return "ForceAtlas2{" +
                "outboundAttCompensation=" + outboundAttCompensation +
                ", edgeWeightInfluence=" + edgeWeightInfluence +
                ", jitterTolerance=" + jitterTolerance +
                ", scalingRatio=" + scalingRatio +
                ", gravity=" + gravity +
                ", speed=" + speed +
                ", speedEfficiency=" + speedEfficiency +
                ", outboundAttractionDistribution=" + outboundAttractionDistribution +
                ", adjustSizes=" + adjustSizes +
                ", barnesHutOptimize=" + barnesHutOptimize +
                ", barnesHutTheta=" + barnesHutTheta +
                ", linLogMode=" + linLogMode +
                ", strongGravityMode=" + strongGravityMode +
                ", threadCount=" + threadCount +
                ", currentThreadCount=" + currentThreadCount +
                ", updateBarnesHutIter=" + updateBarnesHutIter +
                ", updateCenter=" + updateCenter +
                ", barnesHutSplits=" + barnesHutSplits +
                '}';
    }

    private static class SpeedTask implements Callable<Double[]> {
        private final Collection<Node> nodes;

        private SpeedTask(final Collection<Node> nodes) {
            this.nodes = nodes;
        }

        @Override
        public Double[] call() {
            double totalSwinging = 0d;  // How much irregular movement
            double totalEffectiveTraction = 0d;  // Hom much useful movement
            for (Node n : nodes) {
                ForceAtlas2LayoutData nLayout = n.getLayoutData();
                if (!n.isFixed()) {
                    double swinging = Math.sqrt(Math.pow(nLayout.getOld_dx() - nLayout.getDx(), 2) + Math.pow(nLayout.getOld_dy() - nLayout.getDy(), 2) + Math.pow(nLayout.getOld_dz() - nLayout.getDz(), 2));
                    totalSwinging += nLayout.getMass() * swinging;   // If the node has a burst change of direction, then it's not converging.
                    totalEffectiveTraction += nLayout.getMass() * 0.5 * Math.sqrt(Math.pow(nLayout.getOld_dx() + nLayout.getDx(), 2) + Math.pow(nLayout.getOld_dy() + nLayout.getDy(), 2)
                            + Math.pow(nLayout.getOld_dz() + nLayout.getDz(), 2));
                }
            }
            return new Double[]{totalSwinging, totalEffectiveTraction};
        }
    }

    private static class InitLayoutTask implements Runnable {
        private final Collection<Node> nodes;

        public InitLayoutTask(final Collection<Node> nodes, Graph graph) {
            this.nodes = nodes;
        }

        @Override
        public void run() {
            for (Node n : nodes) {
                ForceAtlas2LayoutData nLayout = n.getLayoutData();
//                nLayout.mass = 1 + graph.getDegree(n); mass is constant
                nLayout.setOld_dx(nLayout.getDx());
                nLayout.setOld_dy(nLayout.getDy());
                nLayout.setOld_dz(nLayout.getDz());
                nLayout.setDx(0);
                nLayout.setDy(0);
                nLayout.setDz(0);
            }
        }
    }

    private static class ApplyForcesTask implements Callable<Double> {

        private final Collection<Node> nodes;
        private final boolean adjustSizes;
        private final double speed;
        private final boolean useAltSpeed;

        private ApplyForcesTask(final List<Node> nodes, final boolean adjustSizes, final double speed, final boolean useAltSpeed) {
            this.nodes = nodes;
            this.adjustSizes = adjustSizes;
            this.speed = speed;
            this.useAltSpeed = useAltSpeed;
        }

        @Override
        public Double call() {
            double distance = 0;
            if (adjustSizes) {
                // If nodes overlap prevention is active, it's not possible to trust the swinging measure.
                for (Node n : nodes) {
                    ForceAtlas2LayoutData nLayout = n.getLayoutData();
                    if (!n.isFixed()) {

                        // Adaptive auto-speed: the speed of each node is lowered
                        // when the node swings.
                        double swinging = nLayout.getMass() * Math.sqrt((nLayout.getOld_dx() - nLayout.getDx()) * (nLayout.getOld_dx() - nLayout.getDx()) + (nLayout.getOld_dy() - nLayout.getDy()) * (nLayout.getOld_dy() - nLayout.getDy())
                                + (nLayout.getOld_dz() - nLayout.getDz()) * (nLayout.getOld_dz() - nLayout.getDz()));
                        double factor = 0.1 * speed / (1f + Math.sqrt(speed * swinging));

                        double df = Math.sqrt(Math.pow(nLayout.getDx(), 2) + Math.pow(nLayout.getDy(), 2) + Math.pow(nLayout.getDz(), 2));
                        factor = Math.min(factor * df, 10.) / df;

                        double x = n.x() + nLayout.getDx() * factor;
                        double y = n.y() + nLayout.getDy() * factor;
                        double z = n.z() + nLayout.getDz() * factor;
                        distance += Math.sqrt(Math.pow(n.x() - x, 2) + Math.pow(n.y() - y, 2) + Math.pow(n.z() - z, 2));

                        n.setX((float) x);
                        n.setY((float) y);
                        n.setZ((float) z);
                    }
                }
            } else {
                for (Node n : nodes) {
                    ForceAtlas2LayoutData nLayout = n.getLayoutData();
                    if (!n.isFixed()) {

                        // Adaptive auto-speed: the speed of each node is lowered
                        // when the node swings.
                        double swinging = nLayout.getMass() * Math.sqrt((nLayout.getOld_dx() - nLayout.getDx()) * (nLayout.getOld_dx() - nLayout.getDx()) + (nLayout.getOld_dy() - nLayout.getDy()) * (nLayout.getOld_dy() - nLayout.getDy())
                                + (nLayout.getOld_dz() - nLayout.getDz()) * (nLayout.getOld_dz() - nLayout.getDz()));
                        double factor = speed / (1f + Math.sqrt(speed * swinging));

                        if (useAltSpeed) {
                            factor = Math.min(0.1 * factor, 10.0 / Math.sqrt(Math.pow(nLayout.getDx(), 2) + Math.pow(nLayout.getDy(), 2) + Math.pow(nLayout.getDz(), 2)));
                        }

                        double x = n.x() + nLayout.getDx() * factor;
                        double y = n.y() + nLayout.getDy() * factor;
                        double z = n.z() + nLayout.getDz() * factor;
                        distance += Math.sqrt(Math.pow(n.x() - x, 2) + Math.pow(n.y() - y, 2) + Math.pow(n.z() - z, 2));

                        n.setX((float) x);
                        n.setY((float) y);
                        n.setZ((float) z);
                    }
                }
            }
            return distance;
        }
    }

    private static class BarnesHutUpdateCenterTask implements Runnable {
        private final Collection<Region> regions;
        private final boolean recursive;

        private BarnesHutUpdateCenterTask(final Collection<Region> regions, final boolean recursive) {
            this.regions = regions;
            this.recursive = recursive;
        }

        @Override
        public void run() {
            for (Region r : regions) {
                r.updateAllMassAndGeometry(recursive);
            }
        }
    }

    private static class BarnesHutBuildSubRegionTask implements Runnable {
        private final Collection<Region> regions;
        private final boolean recursive;

        private BarnesHutBuildSubRegionTask(final Collection<Region> regions, final boolean recursive) {
            this.regions = regions;
            this.recursive = recursive;
        }

        @Override
        public void run() {
            for (Region r : regions) {
                r.buildSubRegions(recursive);
            }
        }
    }


    private static class AttractionTask implements Runnable {
        private final Collection<Edge> edges;
        private final ForceFactory.AttractionForce Attraction;
        private final boolean isDynamicWeight;
        private final Interval interval;
        private final double edgeWeightInfluence;

        private AttractionTask(
                final Collection<Edge> edges,
                final ForceFactory.AttractionForce attraction,
                final boolean isDynamicWeight,
                final Interval interval,
                final double edgeWeightInfluence)
        {
            this.edges = edges;
            Attraction = attraction;
            this.isDynamicWeight = isDynamicWeight;
            this.interval = interval;
            this.edgeWeightInfluence = edgeWeightInfluence;
        }

        @Override
        public void run() {
            if (edgeWeightInfluence == 0) {
                for (Edge e : edges) {
                    Attraction.apply(e.getSource(), e.getTarget(), 1);
                }
            } else if (edgeWeightInfluence == 1) {
                for (Edge e : edges) {
                    Attraction.apply(e.getSource(), e.getTarget(), getEdgeWeight(e, isDynamicWeight, interval));
                }
            } else {
                for (Edge e : edges) {
                    Attraction.apply(e.getSource(), e.getTarget(), Math.pow(getEdgeWeight(e, isDynamicWeight, interval), edgeWeightInfluence));
                }
            }
        }
    }
}

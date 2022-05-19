/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis.h>

using namespace gmx;

class Mindistances : public TrajectoryAnalysisModule
{
    public:
        Mindistances();

        virtual void initOptions(IOptionsContainer                    *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        class ModuleData;

        std::string                      fnDist_;
        double                           cutoff_;
        int                              ndist_;
        Selection                        refsel_;
        Selection                        sel_;

        AnalysisNeighborhood             nb_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;
};


Mindistances::Mindistances()
    : cutoff_(0.0)
{
    registerAnalysisDataset(&data_, "avedist");
}


void
Mindistances::initOptions(IOptionsContainer                    *options,
                              TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "This is a template for writing your own analysis tools for",
        "Gromacs. The advantage of using Gromacs for this is that you",
        "have access to all information in the topology, and your",
        "program will be able to handle all types of coordinates and",
        "trajectory files supported by Gromacs. In addition,",
        "you get a lot of functionality for free from the trajectory",
        "analysis library, including support for flexible dynamic",
        "selections. Go ahead an try it![PAR]",
        "To get started with implementing your own analysis program,",
        "follow the instructions in the README file provided.",
        "This template implements a simple analysis programs that calculates",
        "average distances from a reference group to one or more",
        "analysis groups."
    };

    settings->setHelpText(desc);
    //options->setDescription(desc);

    options->addOption(FileNameOption("o")
                           .filetype(eftPlot).outputFile()
                           .store(&fnDist_).defaultBasename("avedist")
                           .description("Average distances from reference group"));

    options->addOption(SelectionOption("reference")
                           .store(&refsel_).required()
                           .description("Reference group to calculate distances from"));
    options->addOption(SelectionOption("select")
                           .store(&sel_).required()
                           .description("Group to calculate distances to"));
    options->addOption(DoubleOption("cutoff").store(&cutoff_)
                           .description("Cutoff for distance calculation (0 = no cutoff)"));
    options->addOption(IntegerOption("ndist").store(&ndist_).defaultValue(2)
		       .description("Number of shortest distances to print"));
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}

void
Mindistances::initAnalysis(const TrajectoryAnalysisSettings &settings,
                               const TopologyInformation         & /*top*/)
{
    nb_.setCutoff(cutoff_);

    data_.setColumnCount(0, ndist_);

    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);

    if (!fnDist_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnDist_);
        plotm->setTitle("Shortest distances");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        plotm->setYFormat(10,5);
        data_.addModule(plotm);
    }
}


void
Mindistances::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle         dh     = pdata->dataHandle(data_);
    const Selection           &refsel = pdata->parallelSelection(refsel_);
    const Selection           &sel   = pdata->parallelSelection(sel_);
    const real INF = 1000.0;
    real dist[ndist_];
    int k,m;

    for(k=0;k<ndist_;k++) dist[k]=INF;
    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, refsel);
    dh.startFrame(frnr, fr.time);
    int              nr    = sel.posCount();
    for (int i = 0; i < nr; ++i) {
      SelectionPosition p = sel.position(i);
      real d = nbsearch.minimumDistance(p.x());
      for(k=0;k<ndist_;k++) if(d<dist[k]) break;
      for(m=ndist_-1; m>k; m--) dist[m]=dist[m-1];
      dist[k]=d;
    }
    for(k=0;k<ndist_;k++) dh.setPoint(k, dist[k]);
    dh.finishFrame();
}


void
Mindistances::finishAnalysis(int /*nframes*/)
{
}


void
Mindistances::writeOutput()
{
    // We print out the average of the mean distances for each group.
    for (int g = 0; g < ndist_; ++g)
    {
        fprintf(stderr, "Average mindistance %d: %.3f nm\n",
                g+1, avem_->average(0, g));
    }
}

/*! \brief
 * The main function for the analysis template.
 */
int
main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<Mindistances>(argc, argv);
}

#include "FilterTimeHits.h"
#include <iostream>
#include <cmath>
#include <set>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

#include <marlin/AIDAProcessor.h>

using namespace lcio;
using namespace marlin;

FilterTimeHits aFilterTimeHits;

FilterTimeHits::FilterTimeHits() : Processor("FilterTimeHits")
{
    // --- Processor description:

    _description = "FilterTimeHits selects tracker hits based on their time corrected for the time of flight";

    // --- Processor parameters:

    registerProcessorParameter("TrackerHitInputCollections",
                               "Name of the tracker hit input collections",
                               m_inputTrackerHitsCollNames,
                               {});

    registerProcessorParameter("TrackerSimHitInputCollections",
                               "Name of the tracker simhit input collections",
                               m_inputTrackerSimHitsCollNames,
                               {});

    registerProcessorParameter("TrackerHitInputRelations",
                               "Name of the tracker hit relation collections",
                               m_inputTrackerHitRelNames,
                               {});

    registerProcessorParameter("TrackerHitOutputCollections",
                               "Name of the tracker hit output collections",
                               m_outputTrackerHitsCollNames,
                               {});

    registerProcessorParameter("TrackerSimHitOutputCollections",
                               "Name of the tracker simhit output collections",
                               m_outputTrackerSimHitsCollNames,
                               {});

    registerProcessorParameter("TrackerHitOutputRelations",
                               "Name of the tracker hit relation collections",
                               m_outputTrackerHitRelNames,
                               {});

    registerProcessorParameter("TargetBeta",
                               "Target beta=v/c for hit time of flight correction",
                               m_beta,
                               double(1.0));

    registerProcessorParameter("TimeLowerLimit",
                               "Lower limit on the corrected hit time in ns",
                               m_time_min,
                               double(-90.0));

    registerProcessorParameter("TimeUpperLimit",
                               "Upper limit on the corrected hit time in ns",
                               m_time_max,
                               double(90.0));

    registerProcessorParameter("FillHistograms",
                               "Flag to fill the diagnostic histograms",
                               m_fillHistos,
                               false);
}

void FilterTimeHits::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // --- Print the processor parameters:

    printParameters();

    // --- Initialize the run and event counters:

    _nRun = 0;
    _nEvt = 0;

    // --- Initialize the AIDAProcessor and book the diagnostic histograms:

    AIDAProcessor::histogramFactory(this);

    m_corrected_time = new TH1F("m_corrected_time", "Corrected time of the hit [ps]", 1000, -250., 250.);
}

void FilterTimeHits::processRunHeader(LCRunHeader *)
{
    _nRun++;
}

void FilterTimeHits::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
                         << "   in run:  " << evt->getRunNumber() << std::endl;

    // --- Check whether the number of input and output collections match

    if (m_inputTrackerHitsCollNames.size() != m_inputTrackerSimHitsCollNames.size() ||
        m_inputTrackerHitsCollNames.size() != m_inputTrackerHitRelNames.size())
    {

        std::stringstream err_msg;
        err_msg << "Mismatch between the reco and sim hits input collections"
                << std::endl;

        throw EVENT::Exception(err_msg.str());
    }

    if (m_outputTrackerHitsCollNames.size() != m_outputTrackerSimHitsCollNames.size() ||
        m_outputTrackerHitsCollNames.size() != m_outputTrackerHitRelNames.size())
    {

        std::stringstream err_msg;
        err_msg << "Mismatch between the reco and sim hits output collections"
                << std::endl;

        throw EVENT::Exception(err_msg.str());
    }

    streamlog_out(DEBUG) << "   passed container size checks" << std::endl;

    // --- Get the input hit collections and create the corresponding output collections:

    const unsigned int nTrackerHitCol = m_inputTrackerHitsCollNames.size();
    std::vector<LCCollection *> inputHitColls(nTrackerHitCol);
    std::vector<LCCollection *> inputSimHitColls(nTrackerHitCol);
    std::vector<LCCollection *> inputHitRels(nTrackerHitCol);

    std::vector<LCCollectionVec *> outputTrackerHitColls(nTrackerHitCol);
    std::vector<LCCollectionVec *> outputTrackerSimHitColls(nTrackerHitCol);
    std::vector<LCCollectionVec *> outputTrackerHitRels(nTrackerHitCol);

    for (unsigned int icol = 0; icol < nTrackerHitCol; ++icol)
    {

        // get the reco hits
        try
        {
            inputHitColls[icol] = evt->getCollection(m_inputTrackerHitsCollNames[icol]);
        }
        catch (lcio::DataNotAvailableException &e)
        {
            streamlog_out(WARNING) << m_inputTrackerHitsCollNames[icol]
                                   << " collection not available" << std::endl;
            continue;
        }

        // get the sim hits
        try
        {
            inputSimHitColls[icol] = evt->getCollection(m_inputTrackerSimHitsCollNames[icol]);
        }
        catch (lcio::DataNotAvailableException &e)
        {
            streamlog_out(WARNING) << m_inputTrackerSimHitsCollNames[icol]
                                   << " collection not available" << std::endl;
            continue;
        }

        // get the reco-sim relations
        try
        {
            inputHitRels[icol] = evt->getCollection(m_inputTrackerHitRelNames[icol]);
        }
        catch (lcio::DataNotAvailableException &e)
        {
            streamlog_out(WARNING) << m_inputTrackerHitRelNames[icol]
                                   << " collection not available" << std::endl;
            continue;
        }

        // reco hit output collections
        std::string encoderString = inputHitColls[icol]->getParameters().getStringVal("CellIDEncoding");
        outputTrackerHitColls[icol] = new LCCollectionVec(inputHitColls[icol]->getTypeName());
        outputTrackerHitColls[icol]->parameters().setValue("CellIDEncoding", encoderString);
        LCFlagImpl lcFlag(inputHitColls[icol]->getFlag());
        outputTrackerHitColls[icol]->setFlag(lcFlag.getFlag());

        // sim hit output collections
        outputTrackerSimHitColls[icol] = new LCCollectionVec(inputSimHitColls[icol]->getTypeName());
        outputTrackerSimHitColls[icol]->parameters().setValue("CellIDEncoding", encoderString);
        LCFlagImpl lcFlag_sim(inputSimHitColls[icol]->getFlag());
        outputTrackerSimHitColls[icol]->setFlag(lcFlag_sim.getFlag());

        // reco-sim relation output collections
        outputTrackerHitRels[icol] = new LCCollectionVec(inputHitRels[icol]->getTypeName());
        LCFlagImpl lcFlag_rel(inputHitRels[icol]->getFlag());
        outputTrackerHitRels[icol]->setFlag(lcFlag_rel.getFlag());
    }

    // --- Loop over the tracker hits and select hits inside the chosen time window:
    std::vector<std::set<int>> hits_to_save(nTrackerHitCol);

    for (unsigned int icol = 0; icol < inputHitColls.size(); ++icol)
    {

        LCCollection *hit_col = inputHitColls[icol];
        if (!hit_col)
            continue;

        for (int ihit = 0; ihit < hit_col->getNumberOfElements(); ++ihit)
        {

            if (TrackerHitPlane *hit = dynamic_cast<TrackerHitPlane *>(hit_col->getElementAt(ihit))){
                // Skipping the hit if its time is outside the acceptance time window
                double hitT = hit->getTime();

                dd4hep::rec::Vector3D pos = hit->getPosition();
                double hitR = pos.r();

                // Correcting for the propagation time
                double dt = hitR / (TMath::C() * m_beta / 1e6);
                hitT -= dt;
                streamlog_out(DEBUG3) << "corrected hit at R: " << hitR << " mm by propagation time: " << dt << " ns to T: " << hitT << " ns" << std::endl;

                //Apply time window selection
                if (hitT < m_time_min || hitT > m_time_max)
                {
                    streamlog_out(DEBUG4) << "hit at T: " << hitT << " ns is rejected by timing cuts" << std::endl;
                    continue;
                }

                hits_to_save[icol].insert(ihit);

                if (m_fillHistos)
                    m_corrected_time->Fill(hitT);
            }

        } // ihit loop

    } // icol loop

    // --- Add the filtered hits to the output collections:

    for (unsigned int icol = 0; icol < inputHitColls.size(); ++icol)
    {

        for (auto &ihit : hits_to_save[icol])
        {

            TrackerHitPlane *hit = dynamic_cast<TrackerHitPlane *>(inputHitColls[icol]->getElementAt(ihit));
            TrackerHitPlaneImpl *hit_new = new TrackerHitPlaneImpl();

            hit_new->setCellID0(hit->getCellID0());
            hit_new->setCellID1(hit->getCellID1());
            hit_new->setType(hit->getType());
            hit_new->setPosition(hit->getPosition());
            hit_new->setU(hit->getU());
            hit_new->setV(hit->getV());
            hit_new->setdU(hit->getdU());
            hit_new->setdV(hit->getdV());
            hit_new->setEDep(hit->getEDep());
            hit_new->setEDepError(hit->getEDepError());
            hit_new->setTime(hit->getTime());
            hit_new->setQuality(hit->getQuality());

	    const lcio::LCObjectVec &rawHits = hit->getRawHits();
	    for (size_t j=0; j<rawHits.size(); ++j) {
	      lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit*>( rawHits[j] );
	      hit_new->rawHits().push_back(hitConstituent);
	    }
            outputTrackerHitColls[icol]->addElement(hit_new);

            LCRelation *rel = dynamic_cast<LCRelation *>(inputHitRels[icol]->getElementAt(ihit));

            SimTrackerHit *simhit = dynamic_cast<SimTrackerHit *>(rel->getTo());
            SimTrackerHitImpl *simhit_new = new SimTrackerHitImpl();

            simhit_new->setCellID0(simhit->getCellID0());
            simhit_new->setCellID1(simhit->getCellID1());
            simhit_new->setPosition(simhit->getPosition());
            simhit_new->setEDep(simhit->getEDep());
            simhit_new->setTime(simhit->getTime());
            simhit_new->setMCParticle(simhit->getMCParticle());
            simhit_new->setMomentum(simhit->getMomentum());
            simhit_new->setPathLength(simhit->getPathLength());
            simhit_new->setQuality(simhit->getQuality());
            simhit_new->setOverlay(simhit->isOverlay());
            simhit_new->setProducedBySecondary(simhit->isProducedBySecondary());

            outputTrackerSimHitColls[icol]->addElement(simhit_new);

            LCRelationImpl *rel_new = new LCRelationImpl();

            rel_new->setFrom(hit_new);
            rel_new->setTo(simhit_new);
            rel_new->setWeight(rel->getWeight());

            outputTrackerHitRels[icol]->addElement(rel_new);

        } // ihit loop

        streamlog_out(MESSAGE) << " " << hits_to_save[icol].size() << " hits added to the collections: "
                               << m_outputTrackerHitsCollNames[icol] << ", "
                               << m_outputTrackerSimHitsCollNames[icol] << ", "
                               << m_outputTrackerHitRelNames[icol] << std::endl;

        evt->addCollection(outputTrackerHitColls[icol], m_outputTrackerHitsCollNames[icol]);
        evt->addCollection(outputTrackerSimHitColls[icol], m_outputTrackerSimHitsCollNames[icol]);
        evt->addCollection(outputTrackerHitRels[icol], m_outputTrackerHitRelNames[icol]);

        streamlog_out(DEBUG5) << " output collection " << m_outputTrackerHitsCollNames[icol] << " of type "
                              << outputTrackerHitColls[icol]->getTypeName() << " added to the event \n"
                              << " output collection " << m_outputTrackerSimHitsCollNames[icol] << " of type "
                              << outputTrackerSimHitColls[icol]->getTypeName() << " added to the event \n"
                              << " output collection " << m_outputTrackerHitRelNames[icol] << " of type "
                              << outputTrackerHitRels[icol]->getTypeName() << " added to the event  "
                              << std::endl;

    } // icol loop

    streamlog_out(DEBUG) << "   Event processed " << std::endl;

    _nEvt++;
}

void FilterTimeHits::check(LCEvent *)
{
}

void FilterTimeHits::end()
{

    std::cout << "FilterTimeHits::end()  " << name()
              << " processed " << _nEvt << " events in " << _nRun << " runs "
              << std::endl;
}

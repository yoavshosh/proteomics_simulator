# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 10:34:57 2019

@author: shosh

This script prepare the mqpar.xml file required for MaxQuant command and run the MQ cmd for analysis as well.
It is according to the template of the 1.6.10.43 version 
"""
import argparse
import subprocess
import time
import os
import shutil

stages_list = ['Configuring',
               'Testing fasta files',
               'Testing raw files',
               'Feature detection',
               'Calculating peak properties',
               'Combining apl files for first search',
               'Preparing searches',
               'MS/MS first search',
               'Read search results for recalibration',
               'Mass recalibration',
               'MS/MS preparation for main search',
               'Combining apl files for main search',
               'MS/MS main search',
               'Preparing combined folder',
               'Calculating masses',
               'Correcting errors',
               'Reading search engine results',
               'Preparing reverse hits', #this is an itraq stage
               'Finish search engine results',
               'Filter identifications (MS/MS)',
               'Applying FDR',
               'Assembling second peptide MS/MS',
               'Combining second peptide files',
               'Second peptide search',
               'Reading search engine results (SP)',
               'Finish search engine results (SP)',
               'Filtering identifications (SP)',
               'Apply FDR (SP)',
               'Re-quantification',
               'Reporter quantification',
               'Prepare protein assembly',
               'Assembling proteins',
               'Assembling unidentified peptides',
               'Finish protein assembly',
               'Updating identifications',
               'Label-free normalization',  #this is an lfq stage
               'Label-free quantification', #this is an lfq stage
               'Label-free collect', #this is an lfq stage
               'Estimating complexity',
               'Prepare writing tables',
               'Writing tables',
               'Finish writing tables']



def write_mqpar_file(raw_files_path,raw_files,fasta_files,fdr,min_pep_len,max_pep_mass,max_mc,enzymes,quantification,phospho,mqpar_name):
    """
    This function writes the mapar file requiered by MaxQuant
    raw_files_path - is the directory of all rawfiles where mqpar is going to be created
    raw_files and fasta_files are list of full paths to raw files and fasta file respectively
    other parameters are just some of the parameters we can set on maxquant
    """
    n_raw_files = len(raw_files)
    
    with open(raw_files_path + 'mqpar_' + mqpar_name + '.xml', 'w') as f:

        f.write('<?xml version="1.0" encoding="utf-8"?>\n')
        f.write('<MaxQuantParams xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n')
        
        f.write('   <fastaFiles>\n')
        for file in fasta_files:
            f.write('      <FastaFileInfo>\n')
            f.write('         <fastaFilePath>'+file+'</fastaFilePath>\n')
            f.write('         <identifierParseRule>>([^\s]*)</identifierParseRule>\n')
            f.write('         <descriptionParseRule>>(.*)</descriptionParseRule>\n')
            f.write('         <taxonomyParseRule></taxonomyParseRule>\n')
            f.write('         <variationParseRule></variationParseRule>\n')
            f.write('         <modificationParseRule></modificationParseRule>\n')
            f.write('         <taxonomyId></taxonomyId>\n')
            f.write('      </FastaFileInfo>\n')
        f.write('   </fastaFiles>\n')
        
        f.write('   <fastaFilesProteogenomics>\n')
        f.write('   </fastaFilesProteogenomics>\n')
        f.write('   <fastaFilesFirstSearch>\n')
        f.write('   </fastaFilesFirstSearch>\n')
        f.write('   <fixedSearchFolder></fixedSearchFolder>\n')
        f.write('   <andromedaCacheSize>350000</andromedaCacheSize>\n')
        f.write('   <advancedRatios>True</advancedRatios>\n')
        f.write('   <pvalThres>0.005</pvalThres>\n')
        f.write('   <neucodeRatioBasedQuantification>False</neucodeRatioBasedQuantification>\n')
        f.write('   <neucodeStabilizeLargeRatios>False</neucodeStabilizeLargeRatios>\n')
        f.write('   <rtShift>False</rtShift>\n')
        f.write('   <separateLfq>False</separateLfq>\n')
        f.write('   <lfqStabilizeLargeRatios>True</lfqStabilizeLargeRatios>\n')
        f.write('   <lfqRequireMsms>True</lfqRequireMsms>\n')
        f.write('   <decoyMode>revert</decoyMode>\n')
        f.write('   <boxCarMode>all</boxCarMode>\n')
        f.write('   <includeContaminants>True</includeContaminants>\n')
        f.write('   <maxPeptideMass>'+str(max_pep_mass)+'</maxPeptideMass>\n')
        f.write('   <epsilonMutationScore>True</epsilonMutationScore>\n')
        f.write('   <mutatedPeptidesSeparately>True</mutatedPeptidesSeparately>\n')
        f.write('   <proteogenomicPeptidesSeparately>True</proteogenomicPeptidesSeparately>\n')
        f.write('   <minDeltaScoreUnmodifiedPeptides>0</minDeltaScoreUnmodifiedPeptides>\n')
        f.write('   <minDeltaScoreModifiedPeptides>6</minDeltaScoreModifiedPeptides>\n')
        f.write('   <minScoreUnmodifiedPeptides>0</minScoreUnmodifiedPeptides>\n')
        f.write('   <minScoreModifiedPeptides>40</minScoreModifiedPeptides>\n')
        f.write('   <secondPeptide>True</secondPeptide>\n')
        f.write('   <matchBetweenRuns>False</matchBetweenRuns>\n')
        f.write('   <matchUnidentifiedFeatures>False</matchUnidentifiedFeatures>\n')
        f.write('   <matchBetweenRunsFdr>False</matchBetweenRunsFdr>\n')
        f.write('   <dependentPeptides>False</dependentPeptides>\n')
        f.write('   <dependentPeptideFdr>0</dependentPeptideFdr>\n')
        f.write('   <dependentPeptideMassBin>0</dependentPeptideMassBin>\n')
        f.write('   <dependentPeptidesBetweenRuns>False</dependentPeptidesBetweenRuns>\n')
        f.write('   <dependentPeptidesWithinExperiment>False</dependentPeptidesWithinExperiment>\n')
        f.write('   <dependentPeptidesWithinParameterGroup>False</dependentPeptidesWithinParameterGroup>\n')
        f.write('   <dependentPeptidesRestrictFractions>False</dependentPeptidesRestrictFractions>\n')
        f.write('   <dependentPeptidesFractionDifference>0</dependentPeptidesFractionDifference>\n')
        f.write('   <msmsConnection>False</msmsConnection>\n')
        f.write('   <ibaq>False</ibaq>\n')
        f.write('   <top3>False</top3>\n')
        f.write('   <independentEnzymes>False</independentEnzymes>\n')
        f.write('   <useDeltaScore>False</useDeltaScore>\n')
        f.write('   <splitProteinGroupsByTaxonomy>False</splitProteinGroupsByTaxonomy>\n')
        f.write('   <taxonomyLevel>Species</taxonomyLevel>\n')
        f.write('   <avalon>False</avalon>\n')
        f.write('   <nModColumns>3</nModColumns>\n')
        f.write('   <ibaqLogFit>False</ibaqLogFit>\n')
        f.write('   <razorProteinFdr>True</razorProteinFdr>\n')
        f.write('   <deNovoSequencing>False</deNovoSequencing>\n')
        f.write('   <deNovoVarMods>True</deNovoVarMods>\n')
        f.write('   <massDifferenceSearch>False</massDifferenceSearch>\n')
        f.write('   <isotopeCalc>False</isotopeCalc>\n')
        f.write('   <writePeptidesForSpectrumFile></writePeptidesForSpectrumFile>\n')
        f.write('   <intensityPredictionsFile>\n')
        f.write('   </intensityPredictionsFile>\n')
        f.write('   <minPepLen>'+str(min_pep_len)+'</minPepLen>\n')
        f.write('   <psmFdrCrosslink>'+str(fdr)+'</psmFdrCrosslink>\n')
        f.write('   <peptideFdr>'+str(fdr)+'</peptideFdr>\n')
        f.write('   <proteinFdr>'+str(fdr)+'</proteinFdr>\n')
        f.write('   <siteFdr>'+str(fdr)+'</siteFdr>\n')
        f.write('   <minPeptideLengthForUnspecificSearch>8</minPeptideLengthForUnspecificSearch>\n')
        f.write('   <maxPeptideLengthForUnspecificSearch>25</maxPeptideLengthForUnspecificSearch>\n')
        f.write('   <useNormRatiosForOccupancy>True</useNormRatiosForOccupancy>\n')
        f.write('   <minPeptides>1</minPeptides>\n')
        f.write('   <minRazorPeptides>1</minRazorPeptides>\n')
        f.write('   <minUniquePeptides>0</minUniquePeptides>\n')
        f.write('   <useCounterparts>False</useCounterparts>\n')
        f.write('   <advancedSiteIntensities>True</advancedSiteIntensities>\n')
        f.write('   <customProteinQuantification>False</customProteinQuantification>\n')
        f.write('   <customProteinQuantificationFile></customProteinQuantificationFile>\n')
        f.write('   <minRatioCount>2</minRatioCount>\n')
        f.write('   <restrictProteinQuantification>True</restrictProteinQuantification>\n')
        f.write('   <restrictMods>\n')
        f.write('      <string>Oxidation (M)</string>\n')
        f.write('      <string>Acetyl (Protein N-term)</string>\n')
        f.write('   </restrictMods>\n')
        f.write('   <matchingTimeWindow>0</matchingTimeWindow>\n')
        f.write('   <matchingIonMobilityWindow>0</matchingIonMobilityWindow>\n')
        f.write('   <alignmentTimeWindow>0</alignmentTimeWindow>\n')
        f.write('   <alignmentIonMobilityWindow>0</alignmentIonMobilityWindow>\n')
        f.write('   <numberOfCandidatesMsms>15</numberOfCandidatesMsms>\n')
        f.write('   <compositionPrediction>0</compositionPrediction>\n')
        f.write('   <quantMode>1</quantMode>\n')
        f.write('   <massDifferenceMods>\n')
        f.write('   </massDifferenceMods>\n')
        f.write('   <mainSearchMaxCombinations>200</mainSearchMaxCombinations>\n')
        f.write('   <writeMsScansTable>False</writeMsScansTable>\n')
        f.write('   <writeMsmsScansTable>True</writeMsmsScansTable>\n')
        f.write('   <writePasefMsmsScansTable>True</writePasefMsmsScansTable>\n')
        f.write('   <writeAccumulatedPasefMsmsScansTable>True</writeAccumulatedPasefMsmsScansTable>\n')
        f.write('   <writeMs3ScansTable>True</writeMs3ScansTable>\n')
        f.write('   <writeAllPeptidesTable>True</writeAllPeptidesTable>\n')
        f.write('   <writeMzRangeTable>True</writeMzRangeTable>\n')
        f.write('   <writeMzTab>False</writeMzTab>\n')
        f.write('   <disableMd5>False</disableMd5>\n')
        f.write('   <cacheBinInds>True</cacheBinInds>\n')
        f.write('   <etdIncludeB>False</etdIncludeB>\n')
        f.write('   <ms2PrecursorShift>0</ms2PrecursorShift>\n')
        f.write('   <complementaryIonPpm>20</complementaryIonPpm>\n')
        f.write('   <variationParseRule></variationParseRule>\n')
        f.write('   <variationMode>none</variationMode>\n')
        f.write('   <useSeriesReporters>False</useSeriesReporters>\n')
        f.write('   <name>session1</name>\n')
        f.write('   <maxQuantVersion>1.6.10.43</maxQuantVersion>\n')
        f.write('   <tempFolder></tempFolder>\n')
        f.write('   <pluginFolder></pluginFolder>\n')
        f.write('   <numThreads>1</numThreads>\n')
        f.write('   <emailAddress></emailAddress>\n')
        f.write('   <smtpHost></smtpHost>\n')
        f.write('   <emailFromAddress></emailFromAddress>\n')
        f.write('   <fixedCombinedFolder></fixedCombinedFolder>\n')
        f.write('   <fullMinMz>-1.79769313486232E+308</fullMinMz>\n')
        f.write('   <fullMaxMz>1.79769313486232E+308</fullMaxMz>\n')
        f.write('   <sendEmail>False</sendEmail>\n')
        f.write('   <ionCountIntensities>False</ionCountIntensities>\n')
        f.write('   <verboseColumnHeaders>False</verboseColumnHeaders>\n')
        f.write('   <calcPeakProperties>False</calcPeakProperties>\n')
        f.write('   <showCentroidMassDifferences>False</showCentroidMassDifferences>\n')
        f.write('   <showIsotopeMassDifferences>False</showIsotopeMassDifferences>\n')
        f.write('   <useDotNetCore>False</useDotNetCore>\n')

        f.write('   <filePaths>\n')
        for file in raw_files:
            f.write('      <string>'+file+'</string>\n')
        f.write('   </filePaths>\n')
        
        f.write('   <experiments>\n')
        for i in range(n_raw_files):
            f.write('      <string></string>\n')      
        f.write('   </experiments>\n')
        
        f.write('   <fractions>\n')
        for i in range(n_raw_files):
            f.write('      <short>32767</short>\n')
        f.write('   </fractions>\n')
           
        f.write('   <ptms>\n')
        for i in range(n_raw_files):
            f.write('      <boolean>False</boolean>\n')
        f.write('   </ptms>\n')
        
        f.write('   <paramGroupIndices>\n')
        for i in range(n_raw_files):
            f.write('      <int>0</int>\n')
        f.write('   </paramGroupIndices>\n')
        
        f.write('   <referenceChannel>\n')
        for i in range(n_raw_files):
            f.write('      <string></string>\n')      
        f.write('   </referenceChannel>\n')  
        
        f.write('   <intensPred>False</intensPred>\n')
        f.write('   <intensPredModelReTrain>False</intensPredModelReTrain>\n')
        f.write('   <parameterGroups>\n')
        f.write('      <parameterGroup>\n')
        f.write('         <msInstrument>0</msInstrument>\n')
        f.write('         <maxCharge>7</maxCharge>\n')
        f.write('         <minPeakLen>2</minPeakLen>\n')
        f.write('         <diaMinPeakLen>2</diaMinPeakLen>\n')
        f.write('         <useMs1Centroids>False</useMs1Centroids>\n')
        f.write('         <useMs2Centroids>False</useMs2Centroids>\n')
        f.write('         <cutPeaks>True</cutPeaks>\n')
        f.write('         <gapScans>1</gapScans>\n')
        f.write('         <minTime>NaN</minTime>\n')
        f.write('         <maxTime>NaN</maxTime>\n')
        f.write('         <matchType>MatchFromAndTo</matchType>\n')
        f.write('         <intensityDetermination>0</intensityDetermination>\n')
        f.write('         <centroidMatchTol>8</centroidMatchTol>\n')
        f.write('         <centroidMatchTolInPpm>True</centroidMatchTolInPpm>\n')
        f.write('         <centroidHalfWidth>35</centroidHalfWidth>\n')
        f.write('         <centroidHalfWidthInPpm>True</centroidHalfWidthInPpm>\n')
        f.write('         <valleyFactor>1.4</valleyFactor>\n')
        f.write('         <isotopeValleyFactor>1.2</isotopeValleyFactor>\n')
        f.write('         <advancedPeakSplitting>False</advancedPeakSplitting>\n')
        f.write('         <intensityThreshold>0</intensityThreshold>\n')
        f.write('         <labelMods>\n')
        f.write('            <string></string>\n')
        f.write('         </labelMods>\n')
        
        if quantification == 'lfq':
            f.write('         <lcmsRunType>Standard</lcmsRunType>\n')
            f.write('         <reQuantify>False</reQuantify>\n')
            f.write('         <lfqMode>1</lfqMode>\n')
        elif quantification in ['itraq','tmt']:
            f.write('         <lcmsRunType>Reporter ion MS2</lcmsRunType>\n')
            f.write('         <reQuantify>False</reQuantify>\n')            
            f.write('         <lfqMode>0</lfqMode>\n')
        
        f.write('         <lfqSkipNorm>False</lfqSkipNorm>\n')
        f.write('         <lfqMinEdgesPerNode>3</lfqMinEdgesPerNode>\n')
        f.write('         <lfqAvEdgesPerNode>6</lfqAvEdgesPerNode>\n')
        f.write('         <lfqMaxFeatures>100000</lfqMaxFeatures>\n')
        f.write('         <neucodeMaxPpm>0</neucodeMaxPpm>\n')
        f.write('         <neucodeResolution>0</neucodeResolution>\n')
        f.write('         <neucodeResolutionInMda>False</neucodeResolutionInMda>\n')
        f.write('         <neucodeInSilicoLowRes>False</neucodeInSilicoLowRes>\n')
        f.write('         <fastLfq>True</fastLfq>\n')
        f.write('         <lfqRestrictFeatures>False</lfqRestrictFeatures>\n')
        f.write('         <lfqMinRatioCount>2</lfqMinRatioCount>\n')
        f.write('         <maxLabeledAa>0</maxLabeledAa>\n')
        f.write('         <maxNmods>5</maxNmods>\n')
        f.write('         <maxMissedCleavages>'+str(max_mc)+'</maxMissedCleavages>\n')
        f.write('         <multiplicity>1</multiplicity>\n')
        f.write('         <enzymeMode>0</enzymeMode>\n')
        f.write('         <complementaryReporterType>0</complementaryReporterType>\n')
        f.write('         <reporterNormalization>0</reporterNormalization>\n')
        f.write('         <neucodeIntensityMode>0</neucodeIntensityMode>\n')
        f.write('         <fixedModifications>\n')
        f.write('            <string>Carbamidomethyl (C)</string>\n')
        f.write('         </fixedModifications>\n')
        
        f.write('         <enzymes>\n')
        for enz in enzymes:
            f.write('            <string>'+enz+'</string>\n')
        f.write('         </enzymes>\n')
        
        f.write('         <enzymesFirstSearch>\n')
        f.write('         </enzymesFirstSearch>\n')
        f.write('         <enzymeModeFirstSearch>0</enzymeModeFirstSearch>\n')
        f.write('         <useEnzymeFirstSearch>False</useEnzymeFirstSearch>\n')
        f.write('         <useVariableModificationsFirstSearch>False</useVariableModificationsFirstSearch>\n')
        f.write('         <variableModifications>\n')
        f.write('            <string>Oxidation (M)</string>\n')
        f.write('            <string>Acetyl (Protein N-term)</string>\n')
        if phospho:
            f.write('            <string>Phospho (STY)</string>\n')
        f.write('         </variableModifications>\n')
        f.write('         <useMultiModification>False</useMultiModification>\n')
        f.write('         <multiModifications>\n')
        f.write('         </multiModifications>\n')
        
        f.write('         <isobaricLabels>\n')
        if quantification == 'itraq':
            for channel_params in [('iTRAQ8plex-Lys113','iTRAQ8plex-Nter113','0','0','0','0','False'),
                                   ('iTRAQ8plex-Lys114','iTRAQ8plex-Nter114','0','0','0','0','False'),
                                   ('iTRAQ8plex-Lys115','iTRAQ8plex-Nter115','0','0','0','0','False'),
                                   ('iTRAQ8plex-Lys116','iTRAQ8plex-Nter116','0','0','0','0','False'),
                                   ('iTRAQ8plex-Lys117','iTRAQ8plex-Nter117','0','0','0','0','False'),
                                   ('iTRAQ8plex-Lys118','iTRAQ8plex-Nter118','0','0','0','0','False'),
                                   ('iTRAQ8plex-Lys119','iTRAQ8plex-Nter119','0','0','0','0','False'),
                                   ('iTRAQ8plex-Lys121','iTRAQ8plex-Nter121','0','0','0','0','False')]:
                f.write('            <IsobaricLabelInfo>\n')
                f.write('               <internalLabel>'+channel_params[0]+'</internalLabel>\n')
                f.write('               <terminalLabel>'+channel_params[1]+'</terminalLabel>\n')
                f.write('               <correctionFactorM2>'+channel_params[2]+'</correctionFactorM2>\n')
                f.write('               <correctionFactorM1>'+channel_params[3]+'</correctionFactorM1>\n')
                f.write('               <correctionFactorP1>'+channel_params[4]+'</correctionFactorP1>\n')
                f.write('               <correctionFactorP2>'+channel_params[5]+'</correctionFactorP2>\n')
                f.write('               <tmtLike>'+channel_params[6]+'</tmtLike>\n')
                f.write('            </IsobaricLabelInfo>\n')
        
        elif quantification == 'tmt':
            for channel_params in [('TMT10plex-Lys126C','TMT10plex-Nter126C','0','0','0','0','True'),
                                   ('TMT10plex-Lys127N','TMT10plex-Nter127N','0','0','0','0','True'),
                                   ('TMT10plex-Lys127C','TMT10plex-Nter127C','0','0','0','0','True'),
                                   ('TMT10plex-Lys128N','TMT10plex-Nter128N','0','0','0','0','True'),
                                   ('TMT10plex-Lys128C','TMT10plex-Nter128C','0','0','0','0','True'),
                                   ('TMT10plex-Lys129N','TMT10plex-Nter129N','0','0','0','0','True'),
                                   ('TMT10plex-Lys129C','TMT10plex-Nter129C','0','0','0','0','True'),
                                   ('TMT10plex-Lys130N','TMT10plex-Nter130N','0','0','0','0','True'),
                                   ('TMT10plex-Lys130C','TMT10plex-Nter130C','0','0','0','0','True'),
                                   ('TMT10plex-Lys131N','TMT10plex-Nter131N','0','0','0','0','True'),
                                   ('TMT11plex-Lys131C','TMT11plex-Nter131C','0','0','0','0','True')]:
                f.write('            <IsobaricLabelInfo>\n')
                f.write('               <internalLabel>'+channel_params[0]+'</internalLabel>\n')
                f.write('               <terminalLabel>'+channel_params[1]+'</terminalLabel>\n')
                f.write('               <correctionFactorM2>'+channel_params[2]+'</correctionFactorM2>\n')
                f.write('               <correctionFactorM1>'+channel_params[3]+'</correctionFactorM1>\n')
                f.write('               <correctionFactorP1>'+channel_params[4]+'</correctionFactorP1>\n')
                f.write('               <correctionFactorP2>'+channel_params[5]+'</correctionFactorP2>\n')
                f.write('               <tmtLike>'+channel_params[6]+'</tmtLike>\n')
                f.write('            </IsobaricLabelInfo>\n')
        f.write('         </isobaricLabels>\n')
        
        f.write('         <neucodeLabels>\n')
        f.write('         </neucodeLabels>\n')
        f.write('         <variableModificationsFirstSearch>\n')
        f.write('         </variableModificationsFirstSearch>\n')
        f.write('         <hasAdditionalVariableModifications>False</hasAdditionalVariableModifications>\n')
        f.write('         <additionalVariableModifications>\n')
        f.write('         </additionalVariableModifications>\n')
        f.write('         <additionalVariableModificationProteins>\n')
        f.write('         </additionalVariableModificationProteins>\n')
        f.write('         <doMassFiltering>True</doMassFiltering>\n')
        f.write('         <firstSearchTol>20</firstSearchTol>\n')
        f.write('         <mainSearchTol>4.5</mainSearchTol>\n')
        f.write('         <searchTolInPpm>True</searchTolInPpm>\n')
        f.write('         <isotopeMatchTol>2</isotopeMatchTol>\n')
        f.write('         <isotopeMatchTolInPpm>True</isotopeMatchTolInPpm>\n')
        f.write('         <isotopeTimeCorrelation>0.6</isotopeTimeCorrelation>\n')
        f.write('         <theorIsotopeCorrelation>0.6</theorIsotopeCorrelation>\n')
        f.write('         <checkMassDeficit>True</checkMassDeficit>\n')
        f.write('         <recalibrationInPpm>True</recalibrationInPpm>\n')
        f.write('         <intensityDependentCalibration>False</intensityDependentCalibration>\n')
        f.write('         <minScoreForCalibration>70</minScoreForCalibration>\n')
        f.write('         <matchLibraryFile>False</matchLibraryFile>\n')
        f.write('         <libraryFile></libraryFile>\n')
        f.write('         <matchLibraryMassTolPpm>0</matchLibraryMassTolPpm>\n')
        f.write('         <matchLibraryTimeTolMin>0</matchLibraryTimeTolMin>\n')
        f.write('         <matchLabelTimeTolMin>0</matchLabelTimeTolMin>\n')
        
        if quantification == 'lfq':
            f.write('         <reporterMassTolerance>NaN</reporterMassTolerance>\n')
            f.write('         <reporterPif>NaN</reporterPif>\n')
            f.write('         <filterPif>False</filterPif>\n')
            f.write('         <reporterFraction>NaN</reporterFraction>\n')
            f.write('         <reporterBasePeakRatio>NaN</reporterBasePeakRatio>\n')
        elif quantification in ['itraq','tmt']:
            f.write('         <reporterMassTolerance>0.003</reporterMassTolerance>\n')
            f.write('         <reporterPif>0</reporterPif>\n')
            f.write('         <filterPif>False</filterPif>\n')
            f.write('         <reporterFraction>0</reporterFraction>\n')
            f.write('         <reporterBasePeakRatio>0</reporterBasePeakRatio>\n')
            
        f.write('         <timsHalfWidth>0</timsHalfWidth>\n')
        f.write('         <timsStep>0</timsStep>\n')
        f.write('         <timsResolution>0</timsResolution>\n')
        f.write('         <timsMinMsmsIntensity>0</timsMinMsmsIntensity>\n')
        f.write('         <timsRemovePrecursor>True</timsRemovePrecursor>\n')
        f.write('         <timsIsobaricLabels>False</timsIsobaricLabels>\n')
        f.write('         <timsCollapseMsms>True</timsCollapseMsms>\n')
        f.write('         <crosslinkSearch>False</crosslinkSearch>\n')
        f.write('         <crossLinker></crossLinker>\n')
        f.write('         <minMatchXl>0</minMatchXl>\n')
        f.write('         <minPairedPepLenXl>6</minPairedPepLenXl>\n')
        f.write('         <crosslinkOnlyIntraProtein>False</crosslinkOnlyIntraProtein>\n')
        f.write('         <crosslinkMaxMonoUnsaturated>0</crosslinkMaxMonoUnsaturated>\n')
        f.write('         <crosslinkMaxMonoSaturated>0</crosslinkMaxMonoSaturated>\n')
        f.write('         <crosslinkMaxDiUnsaturated>0</crosslinkMaxDiUnsaturated>\n')
        f.write('         <crosslinkMaxDiSaturated>0</crosslinkMaxDiSaturated>\n')
        f.write('         <crosslinkModifications>\n')
        f.write('         </crosslinkModifications>\n')
        f.write('         <crosslinkFastaFiles>\n')
        f.write('         </crosslinkFastaFiles>\n')
        f.write('         <crosslinkSites>\n')
        f.write('         </crosslinkSites>\n')
        f.write('         <crosslinkNetworkFiles>\n')
        f.write('         </crosslinkNetworkFiles>\n')
        f.write('         <crosslinkMode></crosslinkMode>\n')
        f.write('         <peakRefinement>False</peakRefinement>\n')
        f.write('         <isobaricSumOverWindow>True</isobaricSumOverWindow>\n')
        f.write('         <isobaricWeightExponent>0.75</isobaricWeightExponent>\n')
        f.write('         <diaLibraryType>0</diaLibraryType>\n')
        f.write('         <diaLibraryPath></diaLibraryPath>\n')
        f.write('         <diaPeptidePaths>\n')
        f.write('         </diaPeptidePaths>\n')
        f.write('         <diaEvidencePaths>\n')
        f.write('         </diaEvidencePaths>\n')
        f.write('         <diaMsmsPaths>\n')
        f.write('         </diaMsmsPaths>\n')
        f.write('         <diaInitialPrecMassTolPpm>20</diaInitialPrecMassTolPpm>\n')
        f.write('         <diaInitialFragMassTolPpm>20</diaInitialFragMassTolPpm>\n')
        f.write('         <diaCorrThresholdFeatureClustering>0.85</diaCorrThresholdFeatureClustering>\n')
        f.write('         <diaPrecTolPpmFeatureClustering>2</diaPrecTolPpmFeatureClustering>\n')
        f.write('         <diaFragTolPpmFeatureClustering>2</diaFragTolPpmFeatureClustering>\n')
        f.write('         <diaScoreN>7</diaScoreN>\n')
        f.write('         <diaMinScore>2.99</diaMinScore>\n')
        f.write('         <diaPrecursorQuant>False</diaPrecursorQuant>\n')
        f.write('         <diaDiaTopNFragmentsForQuant>3</diaDiaTopNFragmentsForQuant>\n')
        f.write('      </parameterGroup>\n')
        f.write('   </parameterGroups>\n')
        f.write('   <msmsParamsArray>\n')
        f.write('      <msmsParams>\n')
        f.write('         <Name>FTMS</Name>\n')
        f.write('         <MatchTolerance>20</MatchTolerance>\n')
        f.write('         <MatchToleranceInPpm>True</MatchToleranceInPpm>\n')
        f.write('         <DeisotopeTolerance>7</DeisotopeTolerance>\n')
        f.write('         <DeisotopeToleranceInPpm>True</DeisotopeToleranceInPpm>\n')
        f.write('         <DeNovoTolerance>10</DeNovoTolerance>\n')
        f.write('         <DeNovoToleranceInPpm>True</DeNovoToleranceInPpm>\n')
        f.write('         <Deisotope>True</Deisotope>\n')
        f.write('         <Topx>12</Topx>\n')
        f.write('         <TopxInterval>100</TopxInterval>\n')
        f.write('         <HigherCharges>True</HigherCharges>\n')
        f.write('         <IncludeWater>True</IncludeWater>\n')
        f.write('         <IncludeAmmonia>True</IncludeAmmonia>\n')
        f.write('         <DependentLosses>True</DependentLosses>\n')
        f.write('         <Recalibration>False</Recalibration>\n')
        f.write('      </msmsParams>\n')
        f.write('      <msmsParams>\n')
        f.write('         <Name>ITMS</Name>\n')
        f.write('         <MatchTolerance>0.5</MatchTolerance>\n')
        f.write('         <MatchToleranceInPpm>False</MatchToleranceInPpm>\n')
        f.write('         <DeisotopeTolerance>0.15</DeisotopeTolerance>\n')
        f.write('         <DeisotopeToleranceInPpm>False</DeisotopeToleranceInPpm>\n')
        f.write('         <DeNovoTolerance>0.25</DeNovoTolerance>\n')
        f.write('         <DeNovoToleranceInPpm>False</DeNovoToleranceInPpm>\n')
        f.write('         <Deisotope>False</Deisotope>\n')
        f.write('         <Topx>8</Topx>\n')
        f.write('         <TopxInterval>100</TopxInterval>\n')
        f.write('         <HigherCharges>True</HigherCharges>\n')
        f.write('         <IncludeWater>True</IncludeWater>\n')
        f.write('         <IncludeAmmonia>True</IncludeAmmonia>\n')
        f.write('         <DependentLosses>True</DependentLosses>\n')
        f.write('         <Recalibration>False</Recalibration>\n')
        f.write('      </msmsParams>\n')
        f.write('      <msmsParams>\n')
        f.write('         <Name>TOF</Name>\n')
        f.write('         <MatchTolerance>40</MatchTolerance>\n')
        f.write('         <MatchToleranceInPpm>True</MatchToleranceInPpm>\n')
        f.write('         <DeisotopeTolerance>0.01</DeisotopeTolerance>\n')
        f.write('         <DeisotopeToleranceInPpm>False</DeisotopeToleranceInPpm>\n')
        f.write('         <DeNovoTolerance>0.02</DeNovoTolerance>\n')
        f.write('         <DeNovoToleranceInPpm>False</DeNovoToleranceInPpm>\n')
        f.write('         <Deisotope>True</Deisotope>\n')
        f.write('         <Topx>10</Topx>\n')
        f.write('         <TopxInterval>100</TopxInterval>\n')
        f.write('         <HigherCharges>True</HigherCharges>\n')
        f.write('         <IncludeWater>True</IncludeWater>\n')
        f.write('         <IncludeAmmonia>True</IncludeAmmonia>\n')
        f.write('         <DependentLosses>True</DependentLosses>\n')
        f.write('         <Recalibration>False</Recalibration>\n')
        f.write('      </msmsParams>\n')
        f.write('      <msmsParams>\n')
        f.write('         <Name>Unknown</Name>\n')
        f.write('         <MatchTolerance>20</MatchTolerance>\n')
        f.write('         <MatchToleranceInPpm>True</MatchToleranceInPpm>\n')
        f.write('         <DeisotopeTolerance>7</DeisotopeTolerance>\n')
        f.write('         <DeisotopeToleranceInPpm>True</DeisotopeToleranceInPpm>\n')
        f.write('         <DeNovoTolerance>10</DeNovoTolerance>\n')
        f.write('         <DeNovoToleranceInPpm>True</DeNovoToleranceInPpm>\n')
        f.write('         <Deisotope>True</Deisotope>\n')
        f.write('         <Topx>12</Topx>\n')
        f.write('         <TopxInterval>100</TopxInterval>\n')
        f.write('         <HigherCharges>True</HigherCharges>\n')
        f.write('         <IncludeWater>True</IncludeWater>\n')
        f.write('         <IncludeAmmonia>True</IncludeAmmonia>\n')
        f.write('         <DependentLosses>True</DependentLosses>\n')
        f.write('         <Recalibration>False</Recalibration>\n')
        f.write('      </msmsParams>\n')
        f.write('   </msmsParamsArray>\n')
        f.write('   <fragmentationParamsArray>\n')
        f.write('      <fragmentationParams>\n')
        f.write('         <Name>CID</Name>\n')
        f.write('         <Connected>False</Connected>\n')
        f.write('         <ConnectedScore0>1</ConnectedScore0>\n')
        f.write('         <ConnectedScore1>1</ConnectedScore1>\n')
        f.write('         <ConnectedScore2>1</ConnectedScore2>\n')
        f.write('         <InternalFragments>False</InternalFragments>\n')
        f.write('         <InternalFragmentWeight>1</InternalFragmentWeight>\n')
        f.write('         <InternalFragmentAas>KRH</InternalFragmentAas>\n')
        f.write('      </fragmentationParams>\n')
        f.write('      <fragmentationParams>\n')
        f.write('         <Name>HCD</Name>\n')
        f.write('         <Connected>False</Connected>\n')
        f.write('         <ConnectedScore0>1</ConnectedScore0>\n')
        f.write('         <ConnectedScore1>1</ConnectedScore1>\n')
        f.write('         <ConnectedScore2>1</ConnectedScore2>\n')
        f.write('         <InternalFragments>False</InternalFragments>\n')
        f.write('         <InternalFragmentWeight>1</InternalFragmentWeight>\n')
        f.write('         <InternalFragmentAas>KRH</InternalFragmentAas>\n')
        f.write('      </fragmentationParams>\n')
        f.write('      <fragmentationParams>\n')
        f.write('         <Name>ETD</Name>\n')
        f.write('         <Connected>False</Connected>\n')
        f.write('         <ConnectedScore0>1</ConnectedScore0>\n')
        f.write('         <ConnectedScore1>1</ConnectedScore1>\n')
        f.write('         <ConnectedScore2>1</ConnectedScore2>\n')
        f.write('         <InternalFragments>False</InternalFragments>\n')
        f.write('         <InternalFragmentWeight>1</InternalFragmentWeight>\n')
        f.write('         <InternalFragmentAas>KRH</InternalFragmentAas>\n')
        f.write('      </fragmentationParams>\n')
        f.write('      <fragmentationParams>\n')
        f.write('         <Name>PQD</Name>\n')
        f.write('         <Connected>False</Connected>\n')
        f.write('         <ConnectedScore0>1</ConnectedScore0>\n')
        f.write('         <ConnectedScore1>1</ConnectedScore1>\n')
        f.write('         <ConnectedScore2>1</ConnectedScore2>\n')
        f.write('         <InternalFragments>False</InternalFragments>\n')
        f.write('         <InternalFragmentWeight>1</InternalFragmentWeight>\n')
        f.write('         <InternalFragmentAas>KRH</InternalFragmentAas>\n')
        f.write('      </fragmentationParams>\n')
        f.write('      <fragmentationParams>\n')
        f.write('         <Name>ETHCD</Name>\n')
        f.write('         <Connected>False</Connected>\n')
        f.write('         <ConnectedScore0>1</ConnectedScore0>\n')
        f.write('         <ConnectedScore1>1</ConnectedScore1>\n')
        f.write('         <ConnectedScore2>1</ConnectedScore2>\n')
        f.write('         <InternalFragments>False</InternalFragments>\n')
        f.write('         <InternalFragmentWeight>1</InternalFragmentWeight>\n')
        f.write('         <InternalFragmentAas>KRH</InternalFragmentAas>\n')
        f.write('      </fragmentationParams>\n')
        f.write('      <fragmentationParams>\n')
        f.write('         <Name>ETCID</Name>\n')
        f.write('         <Connected>False</Connected>\n')
        f.write('         <ConnectedScore0>1</ConnectedScore0>\n')
        f.write('         <ConnectedScore1>1</ConnectedScore1>\n')
        f.write('         <ConnectedScore2>1</ConnectedScore2>\n')
        f.write('         <InternalFragments>False</InternalFragments>\n')
        f.write('         <InternalFragmentWeight>1</InternalFragmentWeight>\n')
        f.write('         <InternalFragmentAas>KRH</InternalFragmentAas>\n')
        f.write('      </fragmentationParams>\n')
        f.write('      <fragmentationParams>\n')
        f.write('         <Name>UVPD</Name>\n')
        f.write('         <Connected>False</Connected>\n')
        f.write('         <ConnectedScore0>1</ConnectedScore0>\n')
        f.write('         <ConnectedScore1>1</ConnectedScore1>\n')
        f.write('         <ConnectedScore2>1</ConnectedScore2>\n')
        f.write('         <InternalFragments>False</InternalFragments>\n')
        f.write('         <InternalFragmentWeight>1</InternalFragmentWeight>\n')
        f.write('         <InternalFragmentAas>KRH</InternalFragmentAas>\n')
        f.write('      </fragmentationParams>\n')
        f.write('      <fragmentationParams>\n')
        f.write('         <Name>Unknown</Name>\n')
        f.write('         <Connected>False</Connected>\n')
        f.write('         <ConnectedScore0>1</ConnectedScore0>\n')
        f.write('         <ConnectedScore1>1</ConnectedScore1>\n')
        f.write('         <ConnectedScore2>1</ConnectedScore2>\n')
        f.write('         <InternalFragments>False</InternalFragments>\n')
        f.write('         <InternalFragmentWeight>1</InternalFragmentWeight>\n')
        f.write('         <InternalFragmentAas>KRH</InternalFragmentAas>\n')
        f.write('      </fragmentationParams>\n')
        f.write('   </fragmentationParamsArray>\n')
        f.write('</MaxQuantParams>\n')

      
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Preparation of parameters file for MaxQuant')
    run_parser = parser.add_argument_group('Run mqpar prep')
    run_parser.add_argument('-raw_path', dest='raw_files_path', action='store', required = True, help='Path to raw input directory')
    run_parser.add_argument('-fasta_files', dest='fasta_files', action='store', nargs='+', required = True, help='list of fasta files with spaces between')
    run_parser.add_argument('-fdr', dest='fdr', action='store', default = '0.01', help='FDR for peptides detection - for now it is the fdr of psmFdrCrosslink\proteins\peptides\sites detection')
    run_parser.add_argument('-min_aa', dest='min_pep_len', action='store', default = '7', help='Minimal number of amino acids per peptide')
    run_parser.add_argument('-max_mass', dest='max_pep_mass', action='store', default = '4600', help='Maximal peptides mass')
    run_parser.add_argument('-max_mc', dest='max_mc', action='store', default = '2', help='maximal number of missed cleavages per peptide')
    run_parser.add_argument('-quantification', dest='quantification', action='store', default = 'lfq', help='quantification method used in MS/MS')
    run_parser.add_argument('-phospho', dest='phospho', action='store', default = 'False', help='raw files are of phosphoproteomics')
    run_parser.add_argument('-enz', dest='enzymes', action='store', nargs='+', default = ['trypsin'], help='enzymes used for cleavage for peptides search. type unspecific for unspecific')
    run_parser.add_argument('-o', dest='output_name', action='store', default = 'combined', help='mqpar and log files name and combined directory (MQ output) name')
    run_parser.add_argument('-run_mq', dest='run_maxquant', action='store', default = 'True', help='run MaxQuant right after mqpar preparation')
    run_parser.add_argument('-delet_intermediate_stages', dest='delet_intermediate_stages', action='store', default = 'True', help='delet the output of intermediate stages')
    run_parser.add_argument('-zip_when_done', dest='zip_raw_files_when_done', action='store', default = 'True', help='zip raw files when MQ is done')
    arguments = parser.parse_args()
    
    raw_files_path = arguments.raw_files_path
    fasta_files = arguments.fasta_files
    fdr = eval(arguments.fdr)
    min_pep_len = eval(arguments.min_pep_len)
    max_pep_mass = eval(arguments.max_pep_mass)
    max_mc = eval(arguments.max_mc)
    output_name = arguments.output_name
    run_maxquant = eval(arguments.run_maxquant)
    enzymes = arguments.enzymes
    quantification = arguments.quantification
    phospho = eval(arguments.phospho)
    delet_intermediate_stages = eval(arguments.delet_intermediate_stages)
    zip_raw_files_when_done = eval(arguments.zip_raw_files_when_done)
    
    #treat complex samples or samples with unknown quantification method as LFQ
    if quantification in ['lfq+','unknown']:
        quantification = 'lfq'
    
    if enzymes == 'unspecific':
        enzymes = []
    if enzymes == ['trypsin']:
        enzymes = ['Trypsin/P']
    
    files_in_raw_dir = [f for f in os.listdir(raw_files_path) if os.path.isfile(os.path.join(raw_files_path, f))] #all files in raw directory
    raw_files = [raw_files_path + f for f in files_in_raw_dir if (f[-4:]=='.raw' or f[-7:]=='.raw.gz' or f[-4:]=='.RAW' or f[-7:]=='.RAW.gz' or f[-5:]=='.wiff'  or f[-8:]=='.wiff.gz' or f[-5:]=='.uimf' or f[-8:]=='.uimf.gz' or f[-6:]=='.mzxml' or f[-9:]=='.mzxml.gz' or f[-6:]=='.mzXML' or f[-9:]=='.mzXML.gz')] #all raw files - full paths
    
    if quantification not in ['lfq','itraq','tmt']:
        print(output_name + ': ' + quantification + 'Quantification method is not supported by this version of script')
    else:
        if run_maxquant:
            #decompress raw files
            decomp=False
            gz_decompress_cmd = 'gzip -d '
            for f in raw_files:
                if f[-3:]=='.gz':
                    decomp=True
                    gz_decompress_cmd += f+' '
            if decomp:
                print('Decompressing raw files')
                try:
                    p1 = subprocess.Popen(gz_decompress_cmd, shell = True, universal_newlines = True)
                    decomp_failed = False
                except subprocess.CalledProcessError as e:
                    err_file = open(raw_files_path + 'err_log.txt','a')
                    print(e.output.decode())
                    err_file.write(e.output.decode()+'\n')
                    err_file.close()  
                    decomp_failed = True
            else:
                decomp_failed = False
                p1 = 1
                
            #wait until decompression is finished
            if decomp and not decomp_failed: 
                while p1.poll() is None: 
                    time.sleep(60*60)
                
            #build mqpar file and run MQ
            if not decomp_failed:
                time.sleep(10)
                raw_files_decompressed = [f[:-3] if f[-3:]=='.gz' else f for f in raw_files] #all decomressed raw files - full paths
                write_mqpar_file(raw_files_path,raw_files_decompressed,fasta_files,fdr,min_pep_len,max_pep_mass,max_mc, enzymes, quantification, phospho, output_name)
                print('Running MQ for '+output_name)
                mq_cmd = 'nohup mono /private/common/Software/MaxQuant/MaxQuant_1.6.10/MaxQuant/bin/MaxQuantCmd.exe ' + raw_files_path+'mqpar_'+output_name+'.xml' + ' > ' + raw_files_path+'nohup_'+output_name+'.out' 
                try:
                    p2 = subprocess.Popen(mq_cmd, shell = True, universal_newlines = True, stdout=subprocess.PIPE)
                    rc2 = p2.poll()
                except subprocess.CalledProcessError as e:
                    err_file = open(raw_files_path + 'err_log.txt','a')
                    print(e.output.decode())
                    err_file.write(e.output.decode()+'\n')
                    err_file.close()             
        
        
    #waiting until MaxQuant finishes - determined by program stdout to nohup log file
    if run_maxquant:
        time.sleep(10)
        while True:
            content = open(raw_files_path+'nohup_'+output_name+'.out', "r").readlines()
            while any(content[-1]==l for l in ['\n','']): #remove uninformative lines from content
                content==content[:-1]
            if content[-1].rstrip() not in stages_list:
                print('MQ failed. stage: ' +content[-1].rstrip())
                time.sleep(10)
                break
            else:
                if content[-1].rstrip()==stages_list[-1]:
                    print('MQ finished')
                    time.sleep(30)
                    break
                else:
                    time.sleep(60*10)
    
    #renaming output directory
    print('Changing combined folder name to '+'combined_'+output_name)
    os.rename(raw_files_path+'combined',raw_files_path+'combined_'+output_name)
    time.sleep(60)
    
    #rm intermediate outputs
    if delet_intermediate_stages:
        print('Deleting intermediate files and folders')
        for f in raw_files_decompressed:
            f_no_suffix = f.replace('.RAW','').replace('.raw','').replace('.wiff','').replace('.uimf','').replace('.mzxml','').replace('.mzXML','')
            os.remove(f_no_suffix+'.index')
            shutil.rmtree(f_no_suffix)
        time.sleep(180*60)
    
    #gzip all raw files
    if zip_raw_files_when_done:
        gzip_compress_cmd = 'gzip '
        for f in raw_files_decompressed:
            gzip_compress_cmd += f+' '
        gzip_compress_cmd += '&'
        print('Compressing raw files')
        try:
            p3 = subprocess.Popen(gzip_compress_cmd, shell = True, universal_newlines = True)
            while p3.poll() is None:
                time.sleep(30)
        except subprocess.CalledProcessError as e:
            err_file = open(raw_files_path + 'err_log.txt','a')
            print(e.output.decode())
            err_file.write(e.output.decode()+'\n')
            err_file.close()
    
    exit(0)
        
    
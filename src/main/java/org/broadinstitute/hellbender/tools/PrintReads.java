package org.broadinstitute.hellbender.tools;

import java.io.IOException;
import java.nio.file.Path;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;

/**
 * Print reads from a SAM/BAM/CRAM file.
 *
 * {@GATK.walkertype ReadWalker}
 */
@CommandLineProgramProperties(
	summary = "Prints reads from the input SAM/BAM/CRAM file to the SAM/BAM/CRAM file.",
    oneLineSummary = "Print reads in the SAM/BAM/CRAM file",
    programGroup = ReadProgramGroup.class
)
@DocumentedFeature
public final class PrintReads extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public String OUTPUT_STR;

    public Path getOutputPath() {
        return IOUtils.getPath(OUTPUT_STR);
    }

    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(getOutputPath(), true);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}

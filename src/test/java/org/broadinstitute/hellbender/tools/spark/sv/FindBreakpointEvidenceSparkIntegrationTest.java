package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Integration test on the SV pipeline as it exists right now [2017-03-06]
 */
public class FindBreakpointEvidenceSparkIntegrationTest extends CommandLineProgramTest {

    private static final String THIS_TEST_FOLDER = getTestDataDir() + "/" + "spark/sv/IntegrationTests";
    private static final String TEST_BAM_LEFT = THIS_TEST_FOLDER + "/" + "assembly13016.bam";
    private static final String TEST_BAM_RIGHT = THIS_TEST_FOLDER + "/" + "assembly13017.bam";
    private static final String KMER_KILL_LIST = THIS_TEST_FOLDER + "/" + "dummy.kill.kmers";
    private static final String ALIGNER_INDEX_IMG = largeFileTestDir + "human_g1k_v37.20.21.fasta.img";
    private static final File reference = new File(b37_reference_20_21);

    private static final class FindBreakpointEvidenceSparkIntegrationTestArgs {
        final String referenceLoc = b37_reference_20_21;
        final String bamLoc;
        final String kmerIgnoreListLoc;
        final String alignerRefIndexImgLoc;
        final String outputDir;

        FindBreakpointEvidenceSparkIntegrationTestArgs(final String bamLoc,
                                                       final String kmerIgnoreListLoc,
                                                       final String alignerRefIndexImgLoc,
                                                       final String outputDir) {
            this.bamLoc = bamLoc;
            this.kmerIgnoreListLoc = kmerIgnoreListLoc;
            this.alignerRefIndexImgLoc = alignerRefIndexImgLoc;
            this.outputDir = outputDir;
        }

        String getCommandLineNoApiKey() {
            return  " -R " + referenceLoc +
                    " -I " + bamLoc +
                    " -O " + outputDir + "/assemblies.sam" +
                    " --alignerIndexImage " + alignerRefIndexImgLoc +
                    " --kmersToIgnore " + kmerIgnoreListLoc +
                    " --breakpointIntervals " + outputDir + "/intervals" +
                    " --fastqDir " + outputDir + "/fastq";
        }

        String getCommandLine() {
            return  getCommandLineNoApiKey() +
                    " --apiKey " + getGCPTestApiKey();
        }
    }

    @DataProvider(name = "findBreakpointEvidenceSparkIntegrationTest")
    public Object[][] createTestData() {

        List<Object[]> tests = new ArrayList<>();
        final File tempDirLeft = BaseTest.createTempDir("forLeft");
        tempDirLeft.deleteOnExit();
        tests.add(new Object[]{new FindBreakpointEvidenceSparkIntegrationTestArgs(TEST_BAM_LEFT, KMER_KILL_LIST, ALIGNER_INDEX_IMG, tempDirLeft.getAbsolutePath())});
        final File tempDirRight = BaseTest.createTempDir("forRight");
        tempDirRight.deleteOnExit();
        tests.add(new Object[]{new FindBreakpointEvidenceSparkIntegrationTestArgs(TEST_BAM_RIGHT, KMER_KILL_LIST, ALIGNER_INDEX_IMG, tempDirRight.getAbsolutePath())});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "findBreakpointEvidenceSparkIntegrationTest", groups = "spark", enabled = false)
    // TODO: 4/5/17 enable after PR#2444 is in
    public void testFindBreakpointRunnableLocal(final FindBreakpointEvidenceSparkIntegrationTestArgs params) throws IOException {

        @SuppressWarnings("unchecked")
        final List<String> expectedFileNames = new ArrayList<>(Collections.EMPTY_LIST);

        new IntegrationTestSpec(
                new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getString(),
                expectedFileNames)
                .executeTest("testFindBreakpointEvidenceSparkRunnableLocal-", this);
    }

    @Test(dataProvider = "findBreakpointEvidenceSparkIntegrationTest", groups = "spark", enabled = false)
    // TODO: 4/5/17 enable after PR#2444 is in
    public void testFindBreakpointRunnableMiniCluster(final FindBreakpointEvidenceSparkIntegrationTestArgs params) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getArgsArray() );
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            int idx = 0;
            // inputs, to be copied to mini-cluster
            idx = argsToBeModified.indexOf("-R");
            Path path = new Path(workingDirectory, "hdfs.fasta");
            cluster.getFileSystem().copyFromLocalFile(new Path(reference.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("-I");
            path = new Path(workingDirectory, "hdfs.bam");
            File file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--kmersToIgnore");
            path = new Path(workingDirectory, "dummy.kill.kmers");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("-O");
            path = new Path(workingDirectory, "assemblies.sam");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--breakpointIntervals");
            path = new Path(workingDirectory, "intervals");
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--fastqDir");
            path = new Path(workingDirectory, "fastq");
            argsToBeModified.set(idx+1, path.toUri().toString());

            @SuppressWarnings("unchecked")
            final List<String> expectedFileNames = new ArrayList<>(Collections.EMPTY_LIST);

            new IntegrationTestSpec(String.join(" ", argsToBeModified), expectedFileNames)
                    .executeTest("testFindBreakpointEvidenceSparkRunnableMiniCluster-", this);
        });
    }
}

package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.NovelAdjacencyReferenceLocations.EndConnectionType.*;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointComplications.DUPSEQ_REPEAT_UNIT_NA_VALUE;

public class NovelAdjacencyReferenceLocationsUnitTest extends BaseTest{

    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!SVCallerTestDataProvider.testDataInitialized) {
            new SVCallerTestDataProvider();
        }
    }

    private static void seeIfItWorks(final NovelAdjacencyReferenceLocations breakpoints, final NovelAdjacencyReferenceLocations.EndConnectionType expectedEndConnectionType,
                                     final SimpleInterval expectedLeftBreakpoint, final SimpleInterval expectedRightBreakpoint,
                                     final SimpleInterval expectedRepeatUnitRefSpan, final String expectedHomology, final String expectedInsertion,
                                     final int expectedRefDupNum, final int expectedCtgDupNum,
                                     final List<String> expectedTaddupCigarStrings) {

        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, expectedLeftBreakpoint);
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, expectedRightBreakpoint);
        Assert.assertEquals(breakpoints.endConnectionType, expectedEndConnectionType);
        Assert.assertEquals(breakpoints.complication.getHomologyForwardStrandRep(), expectedHomology);
        Assert.assertEquals(breakpoints.complication.getInsertedSequenceForwardStrandRep(), expectedInsertion);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), expectedRepeatUnitRefSpan);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), expectedRefDupNum);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), expectedCtgDupNum);
        Assert.assertEquals(breakpoints.complication.getCigarStringsForDupSeqOnCtg(), expectedTaddupCigarStrings);
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for generic functions on the base class
    // -----------------------------------------------------------------------------------------------
    @Test
    public void testEqualsAndHashCode() throws Exception {

        final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations1 = getBreakpoints("1", "contig-1", "foo");

        final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations2 = getBreakpoints("2", "contig-2", "bar");

        Assert.assertTrue(novelAdjacencyReferenceLocations1.equals(novelAdjacencyReferenceLocations2));
        Assert.assertEquals(novelAdjacencyReferenceLocations1.hashCode(), novelAdjacencyReferenceLocations2.hashCode());
    }

    @Test
    void testKryoSerializer() throws IOException {
        // uses inversion subclass for testing
        final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations1 = getBreakpoints("1", "contig-1", "foo");
        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, novelAdjacencyReferenceLocations1);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final NovelAdjacencyReferenceLocations roundTrip = (NovelAdjacencyReferenceLocations)kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, novelAdjacencyReferenceLocations1);
    }

    private static NovelAdjacencyReferenceLocations getBreakpoints(final String assemblyId, final String contigId, final String insertionMapping) throws IOException{
        final AlignmentRegion region1 = new AlignmentRegion(assemblyId, contigId, new SimpleInterval("20", 10000, 10100), TextCigarCodec.decode("100M"), true, 60, 0, 1, 100);
        final AlignmentRegion region2 = new AlignmentRegion(assemblyId, contigId, new SimpleInterval("20", 20100, 20200), TextCigarCodec.decode("100M"), false, 60, 0, 101, 200);
        final ArrayList<String> insertionMappings = new ArrayList<>();
        insertionMappings.add(insertionMapping);
        final ChimericAlignment breakpoint = new ChimericAlignment(region1, region2, SVCallerTestDataProvider.makeDummySequence(100, (byte)'A'), insertionMappings);
        return new NovelAdjacencyReferenceLocations(breakpoint);
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for inversion
    // -----------------------------------------------------------------------------------------------
    @Test
    public void testGetBreakpoints_5to3Inversion_withSimpleInsertion() throws IOException {

        final NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleInversionWithNovelInsertion._3();
        seeIfItWorks(breakpoints, FIVE_TO_FIVE,
                new SimpleInterval("21", 108569294, 108569294), new SimpleInterval("21", 108569364, 108569364),
                DUPSEQ_REPEAT_UNIT_NA_VALUE, "", "T", 0, 0, Collections.emptyList());
    }

    @Test
    public void testGetAssembledBreakpointFromAlignmentRegionsStrangeLeftBreakpoint() throws Exception {

        final NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3();
        seeIfItWorks(breakpoints, THREE_TO_THREE,
                new SimpleInterval(SVCallerTestDataProvider.chrForLongContig1, 20138006, 20138006),
                new SimpleInterval(SVCallerTestDataProvider.chrForLongContig1, 20152650, 20152650),
                DUPSEQ_REPEAT_UNIT_NA_VALUE, "TGAGGTCAGGAGTTCCTGATCCCATCTTTACTAAAAATACAAAACTTACCCAGGGTGGTTGTGCACACTTGTAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATTGCTTGAACAAGGGAGGAAATGGTTGCAGTGAGCCATGATCATGCCACTGAACCCCAGCCTGGGCAAGAGAGTGAGACTGTCTCAAAAAAAAAAAAAACTGTTTAATTTTTATGAATGCAGGTTTTCTGCAAACACTACACATAACTATGCTAATTGTTCTGAAGTAATAAATAGAAAGCAAGGCACAACTACAGACTCCACTGTTCAGTTTATGCACTGAACTGTTCTTGCTTTTGCAGTGTAAGTATTTCTGCCTGCAAATACTGGATAATTACCTTGGATCATCAGATTTCTATCAAAGGAATTTAGTATCTTTTAGTCTTTATCATTTTGTATTGCTAAATTTATCTGTGTGTTAAGCTTCTGTGTGCTCTTAAAATGAGGTTTTATCTAAACAAACCTGTGTCTACTTTAAAAGACTAAACATGAAAAAACTAAACTTTTCAGAACCAAAAACAAAGCAATAAATCTGAAGTACTAGATAGTCTGGAGTGAGATTTATTTAGCTTTTTT",
                "", 0, 0, Collections.emptyList());
    }

    /**
     *  @see SVCallerTestDataProvider#forSimpleInversionWithHomology(ByteArrayOutputStream)
     */
    @Test
    public void testSimpleInversionWithHomologyBreakpointsIdentification_allFourRepresentations() throws IOException {

        final byte[] homology = "ACACA".getBytes();

        // left flanking forward strand
        final NovelAdjacencyReferenceLocations breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand = SVCallerTestDataProvider.forSimpleInversionWithHom_leftPlus._3();
        seeIfItWorks(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand, FIVE_TO_FIVE,
                new SimpleInterval("20", 200, 200), new SimpleInterval("20", 605, 605),
                DUPSEQ_REPEAT_UNIT_NA_VALUE, new String(homology), "", 0, 0, Collections.emptyList());

        // see if reverse strand changes anything
        Assert.assertEquals(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand, SVCallerTestDataProvider.forSimpleInversionWithHom_leftMinus._3());

        // see if right flanking evidence give the same breakpoint location and homology (up to RC)
        // and see if the two strands give the same result
        final NovelAdjacencyReferenceLocations breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand = SVCallerTestDataProvider.forSimpleInversionWithHom_rightPlus._3();
        Assert.assertEquals(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand.leftJustifiedLeftRefLoc,
                breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.leftJustifiedLeftRefLoc);
        Assert.assertEquals(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand.leftJustifiedRightRefLoc,
                breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.leftJustifiedRightRefLoc);
        Assert.assertEquals(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand.complication.getHomologyForwardStrandRep(),
                new String(SVCallerTestDataProvider.getReverseComplimentCopy(breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.complication.getHomologyForwardStrandRep().getBytes())));
        Assert.assertEquals(breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.endConnectionType, THREE_TO_THREE);
        Assert.assertEquals(breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand, SVCallerTestDataProvider.forSimpleInversionWithHom_rightMinus._3());
    }

    @Test
    public void testGetAssembledBreakpointsFromAlignmentRegionsWithOverlappingAlignmentRegion() throws Exception {
        final byte[] contigSequence = "ACTAGAGCATCTACGTGTTCCTGTGGTTTTGGAGCAAGAGTGATTTGAGTTTCAGAGATTTTTACTAATTCTTCTTCCCCTACCAGAAAAAAAGATCTTACCATTTGAGAGTGAGATGTAAACCCAGCCCTGTCTGACCTGAGTCTGTGCCCTAAGCCTATGCTAAGCCAAGCAGTGCCTGGAGCCACCACAGGTCCACACAATTCGTTAACATGATGAAGCAAGGATGGAAATTGGACAAAATAGTGTGCCTACTGAATCTAAGAATGAAAAATGATTGCACTCCTACTCTGAGTGCTTTGGAGCACTGCCCAGTTGGGCAAAGGGTCAGCGCCTGGGCAGAGGTCCCCACAACCTGGCAGGAGTGTGGTCGGCCACCCTATGGGCCTCCATCATGTGCAGTGACAGCGGGGCTGTCATGTCACCGTGTGGGAGGGCTTGCAGGTGAAGTGGTCTGGGAGGGGTCCCCCAGACAAAGCCAAGGTTCTGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCACCTGCAAGCCCTCCCACACGGTGACATGACAGCCTATAATACAGTTCCGCATGGCCACGTCATACAACCCTGTCATATTGGTGAGCAATTGCTGTGTAGCCAAAGACCCCAAAACTCAAACAGCATTTATTATTATTGCCCCCATGTCTGAGAGTCAGATGTGCATTTGCTGATCTCAGCTTGTTTGAGCTGCTGCAGGGTTGGGGCTCTGCTCCAGGCAGGCTTAGCTGTCACCACATGCACACATACATTCTGGGCCTCTGCTGCGCGCGTCACGTTCACTGAAGATCTTGGGATTGGGAGTTAGGGCGGTGGGAGGGCCCAGCAAAGTCACCTGGCGATGGCAGGGACACAGGGAGGAATGTAGAATGGGGCCGATGATGGGACCCACACGTCTGCAAAGCTGCGGTCTCCTTGAGGGGTGGAGACAGCAACAACTCACCGCACGCGGTGCTTCAGTTCACCATCTCCCTGGGACATTAGGGGGCCCCGTGTTATCTCATTTTGCTCTGGTTTGCATTAGTTTTTTATCACTTCGTAGATGAAGCCACTGACACCCAGAGAGGGAAAGTGGCCTGACCAAGGGCCACAGCAGGGGAGCGAAGGAGCCCCACAGTTCGGCAGGAACACAGCCTCTCCCTGGCTTTCAGGTTCACTGACATCTTCTCATGGCCTCTGTAACTCACCAGGCATCAGGGTGTAGTCCTTAGACCAGTGTCCCACAGCTGCCACAGAGTGGGAGCTCACCATCAGTTATAAGTCACTAGAAAGGCTTTTGGACATTATAAGCTACAATGGAAAATAAGTCATCTGTGGATTTTTGTGACAGATTCCAAAAATTTGAATATTTTGTCTACTTAGGTTTTTGGTTAATTTTATCCTCAAAACTGTTCTGCAGTGATTAAGCTGTACAAACTGCATCATGGGCGAATTGGCATATTCAGAAATGACTGATATTCTTGATTTCAGTTTTTTACTTTGTATGTAGCTCCTCAAGGAAAC".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 23102817, 23103304), TextCigarCodec.decode("487M1006S"), true, 60, 1, 1, 487);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 23103196, 23103238), TextCigarCodec.decode("483S42M968S"), false, 60, 2, 484, 525);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 23103633, 23104603), TextCigarCodec.decode("523S970M"), true, 60, 3, 524, 1493);
        final List<AlignmentRegion> alignmentRegionList = Arrays.asList(region1, region2, region3);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = ChimericAlignment.fromSplitAlignments(new Tuple2<>(alignmentRegionList, contigSequence));
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.regionWithLowerCoordOnContig.contigId, "contig-1");
        Assert.assertEquals(chimericAlignment.regionWithLowerCoordOnContig, region1);
        Assert.assertEquals(chimericAlignment.regionWithHigherCoordOnContig, region3);
        Assert.assertEquals(chimericAlignment.insertionMappings.size(), 1);
        final String expectedInsertionMappingsString = String.join(AlignmentRegion.PACKED_STRING_REP_SEPARATOR, "1", "contig-1", "484", "525", "20", "23103196", "-", "483S42M968S", "60", "2");
        Assert.assertEquals(chimericAlignment.insertionMappings.get(0), expectedInsertionMappingsString);//"1_contig-1_484_525_20_23103196_-_483S42M968S_60_2");
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(chimericAlignment);
        Assert.assertTrue(breakpoints.complication.getHomologyForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getInsertedSequenceForwardStrandRep(), "TGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCA");
    }

    // following might be legacy tests that could be removed but needs time to investigate (Dec.13/2016)
    @Test
    public void testStrandedness_Inversion() throws IOException {

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 10000, 10100), TextCigarCodec.decode("100M"), true, 60, 0, 1, 100);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 20100, 20200), TextCigarCodec.decode("100M"), false, 60, 0, 101, 200);
        final ChimericAlignment breakpoint1 = new ChimericAlignment(region1, region2, SVCallerTestDataProvider.makeDummySequence(100, (byte)'A'), new ArrayList<>());

        Assert.assertNotEquals(NovelAdjacencyReferenceLocations.determineEndConnectionType(breakpoint1), NovelAdjacencyReferenceLocations.EndConnectionType.FIVE_TO_THREE);

        final AlignmentRegion region3 = new AlignmentRegion("4", "contig-7", new SimpleInterval("21", 38343346, 38343483), TextCigarCodec.decode("137M141S"), true, 60, 0, 1, 137);
        final AlignmentRegion region4 = new AlignmentRegion("4", "contig-7", new SimpleInterval("20", 38342908, 38343049), TextCigarCodec.decode("137S141M"), false, 60, 0, 138, 278);
        final ChimericAlignment breakpoint2 = new ChimericAlignment(region3, region4, SVCallerTestDataProvider.makeDummySequence(137+141, (byte)'A'), new ArrayList<>());

        Assert.assertEquals(NovelAdjacencyReferenceLocations.determineEndConnectionType(breakpoint2), NovelAdjacencyReferenceLocations.EndConnectionType.FIVE_TO_FIVE);

        final AlignmentRegion region5 = new AlignmentRegion("3", "contig-7", new SimpleInterval("21", 38343346, 38343483), TextCigarCodec.decode("137M141S"), true, 60, 0, 1, 137);
        final AlignmentRegion region6 = new AlignmentRegion("3", "contig-7", new SimpleInterval("21", 38342908, 38343049), TextCigarCodec.decode("137S141M"), false, 60, 0, 138, 278);
        final ChimericAlignment breakpoint3 = new ChimericAlignment(region5, region6, SVCallerTestDataProvider.makeDummySequence(137+141, (byte)'A'), new ArrayList<>());

        Assert.assertEquals(NovelAdjacencyReferenceLocations.determineEndConnectionType(breakpoint3), NovelAdjacencyReferenceLocations.EndConnectionType.FIVE_TO_FIVE);
    }

    @Test
    public void testGetBreakpoints_5to3Inversion_simple() throws IOException {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", new SimpleInterval("20", 101, 200), TextCigarCodec.decode("100M100S"), true, 60, 0, 1, 100);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", new SimpleInterval("20", 501, 600), TextCigarCodec.decode("100S100M"), false, 60, 0, 101, 200);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, SVCallerTestDataProvider.makeDummySequence(200, (byte)'A'), Collections.emptyList());
        final BreakpointComplications complication = new BreakpointComplications("", "", DUPSEQ_REPEAT_UNIT_NA_VALUE, 0, 0, new ArrayList<>(), false);
        final Tuple2<SimpleInterval, SimpleInterval> breakpoints = NovelAdjacencyReferenceLocations.leftJustifyBreakpoints(chimericAlignment, complication);
        Assert.assertEquals(breakpoints._1(), new SimpleInterval("20", 200, 200));
        Assert.assertEquals(breakpoints._2(), new SimpleInterval("20", 600, 600));
    }

    @Test
    public void testGetBreakpoints_5to3Inversion_withSimpleHomology() throws IOException {
        final AlignmentRegion region1 = new AlignmentRegion("1","1", new SimpleInterval("20", 101, 205), TextCigarCodec.decode("105M100S"), true, 60, 0, 1, 105);
        final AlignmentRegion region2 = new AlignmentRegion("1","1", new SimpleInterval("20", 501, 605), TextCigarCodec.decode("105M100S"), false, 60, 0, 96, 200);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, SVCallerTestDataProvider.makeDummySequence(205, (byte)'A'), Collections.emptyList());
        final BreakpointComplications complication = new BreakpointComplications("ACACA", "", DUPSEQ_REPEAT_UNIT_NA_VALUE, 0, 0, new ArrayList<>(), false);
        final Tuple2<SimpleInterval, SimpleInterval> breakpoints = NovelAdjacencyReferenceLocations.leftJustifyBreakpoints(chimericAlignment, complication);
        Assert.assertEquals(breakpoints._1(), new SimpleInterval("20", 200, 200));
        Assert.assertEquals(breakpoints._2(), new SimpleInterval("20", 605, 605));
    }

    @Test
    public void testGetHomology() {
        final byte[] contigSequence = "ATCGATCGAAAAGCTAGCTA".getBytes();

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 1, 12), TextCigarCodec.decode("12M8S"), true, 60, 1, 1, 12);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 101, 112), TextCigarCodec.decode("8H12M"), false, 60, 1, 9, 20);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertEquals(BreakpointComplications.getHomology(region1, region2, contigSequence), "AAAA");

        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 1, 12), TextCigarCodec.decode("8M"), true, 60, 1, 1, 8);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentRegion region4 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 101, 112), TextCigarCodec.decode("8M"), false, 60, 1, 13, 20);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertTrue(BreakpointComplications.getHomology(region3, region4, contigSequence).isEmpty());
    }

    @Test
    public void testGetInsertedSequence() {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("8", 118873207, 118873739), TextCigarCodec.decode("532M87S"), true, 60, 0, 1, 532);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 175705642, 175705671), TextCigarCodec.decode("518S29M72S"), false, 3, 0, 519, 547);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 118875262, 118875338), TextCigarCodec.decode("543S76M"), false, 60, 0, 544, 619);

        Assert.assertTrue(BreakpointComplications.getInsertedSequence(region3, region1, contigSequence).isEmpty());
        Assert.assertEquals(BreakpointComplications.getInsertedSequence(region1, region3, contigSequence), "GAGATAGAGTC");

        Assert.assertTrue(BreakpointComplications.getInsertedSequence(region2, region1, contigSequence).isEmpty() && BreakpointComplications.getInsertedSequence(region1, region2, contigSequence).isEmpty());
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for complication resolving and breakpoint justification with the inferred complications for insertion and deletion
    // -----------------------------------------------------------------------------------------------
    @Test(expectedExceptions = GATKException.class)
    public void testGetBreakpoints_ExpectException() throws IOException {
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100100), TextCigarCodec.decode("100M"), true, 60, 0, 1 ,100);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100101, 100200), TextCigarCodec.decode("100M"), true, 60, 0, 101 ,200);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, SVCallerTestDataProvider.makeDummySequence(200, (byte)'A'), Collections.emptyList());
        new BreakpointComplications(chimericAlignment);
    }

    /**
     * @see SVCallerTestDataProvider#forSimpleDeletion(ByteArrayOutputStream)
     */
    @Test
    public void testGetBreakpoints_simpleDeletion() throws IOException {

        final NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleDeletion_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVCallerTestDataProvider.forSimpleDeletion_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100060, 100060),
                DUPSEQ_REPEAT_UNIT_NA_VALUE, "", "", 0, 0, Collections.emptyList());
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVCallerTestDataProvider#forSimpleInsertion(ByteArrayOutputStream)
     */
    @Test
    public void testGetBreakpoints_simpleInsertion() throws IOException {

        byte[] insertedSeq  = SVCallerTestDataProvider.makeDummySequence(50, (byte)'C');

        final NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleInsertion_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVCallerTestDataProvider.forSimpleInsertion_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100100, 100100), new SimpleInterval("21", 100100, 100100),
                DUPSEQ_REPEAT_UNIT_NA_VALUE, "", new String(insertedSeq), 0, 0, Collections.emptyList());
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVCallerTestDataProvider#forLongRangeSubstitution()
     */
    @Test
    public void testGetBreakpoints_longRangeSubstitution() throws IOException {

        final byte[] substitution = SVCallerTestDataProvider.makeDummySequence(10, (byte)'C');
        final NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forLongRangeSubstitution_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVCallerTestDataProvider.forLongRangeSubstitution_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100060, 100060),
                DUPSEQ_REPEAT_UNIT_NA_VALUE, "", new String(substitution), 0, 0, Collections.emptyList());
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVCallerTestDataProvider#forDeletionWithHomology(ByteArrayOutputStream)
     */
    @Test
    public void testGetBreakpoints_simpleDeletion_withHomology() throws IOException {

        final byte[] homology = "ATCG".getBytes();

        final NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forDeletionWithHomology_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVCallerTestDataProvider.forDeletionWithHomology_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100078, 100078),
                DUPSEQ_REPEAT_UNIT_NA_VALUE, new String(homology), "", 0, 0, Collections.emptyList());
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVCallerTestDataProvider#forSimpleTandemDuplicationContraction()
     */
    @Test
    public void testGetBreakpoints_tandemDuplication_contraction_simple() throws IOException {

        final NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleTanDupContraction_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVCallerTestDataProvider.forSimpleTanDupContraction_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100050, 100050),
                new SimpleInterval("21", 100041, 100050),
                new String(SVCallerTestDataProvider.makeDummySequence(10, (byte)'C')), "",
                2, 1, Collections.emptyList());
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVCallerTestDataProvider#forSimpleTandemDuplicationExpansion(ByteArrayOutputStream)
     */
    @Test
    public void testGetBreakpoints_tandemDuplication_expansion_simple() throws IOException {

        final NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleTanDupExpansion_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVCallerTestDataProvider.forSimpleTanDupExpansion_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100040, 100040),
                new SimpleInterval("21", 100041, 100050),
                "", "",
                1, 2, Arrays.asList("10M", "10M"));
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVCallerTestDataProvider#forSimpleTandemDuplicationExpansionWithNovelInsertion(ByteArrayOutputStream)
     */
    @Test
    public void testGetBreakpoints_tandemDuplication_expansion_andNovelInsertion() throws IOException {

        final byte[] insertedSeq = "CTCTCTCTCT".getBytes();                                                                           //10
        final byte[] dup = "AAAAGTAAATGTTATAAGAAATCTTAAGTATTATTTTCTTATGTTTCTAGCCTAATAAAGTGCTTTTATTAAAGCACTTTATTTAAAGG".getBytes();    //89

        final NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 25297163, 25297163), new SimpleInterval("21", 25297163, 25297163),
                new SimpleInterval("21", 25297164,25297252),
                "", new String(insertedSeq),
                1, 2, Arrays.asList("89M", "89M"));
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVCallerTestDataProvider#forComplexTandemDuplication()
     */
    @Test
    public void testGetBreakpoints_tandemDuplication_complex() throws IOException {

        final String leftRefFlank       = "TGCCAGGTTACATGGCAAAGAGGGTAGATAT";                                                                    // 31
        final String rightRefFlank      = "TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC";                                                            // 39
        final String firstRepeat        = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA";   // 96
        final String secondRepeat       = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA";   // 96
        final String pseudoHomology     = "GGGCAGCTGTGGA";                                                                                      // 13


        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus._3();

        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), pseudoHomology)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 1);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 2);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        breakpoints = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), pseudoHomology)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 1);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 2);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        // second test: contraction from 2 units to 1 unit with pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat+pseudoHomology)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 1);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312705, 312705));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat+pseudoHomology)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 2);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 1);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312705, 312705));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        // third test: contraction from 3 units to 2 units without pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat+secondRepeat)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 3);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 2);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312705, 312705));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        breakpoints = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat+secondRepeat)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 3);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 2);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312705, 312705));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 2);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 3);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 2);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 3);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);
    }

    /**
     * todo : see below
     * Tests a simple inter-chromosomal translocation event, which can not be handled yet [4/4/2017];
     *   enable (and definitely expand test cases) when such variants can be handled
     * @throws Exception
     */
    @Test(enabled = false)
    public void testGetAssembledBreakpointsFromAlignmentRegions() throws Exception {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 118873207, 118873739), TextCigarCodec.decode("532M87S"), true, 60, 0, 1, 532);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 175705642, 175705671), TextCigarCodec.decode("518S29M72S"), false, 3, 0, 519, 547);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 118875262, 118875338), TextCigarCodec.decode("543S76M"), false, 60, 0, 544, 619);
        final List<AlignmentRegion> alignmentRegionList = Arrays.asList(region1, region2, region3);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = ChimericAlignment.fromSplitAlignments(new Tuple2<>(alignmentRegionList, contigSequence));
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.regionWithLowerCoordOnContig.contigId, "contig-1");
        Assert.assertEquals(chimericAlignment.regionWithLowerCoordOnContig, region1);
        Assert.assertEquals(chimericAlignment.regionWithHigherCoordOnContig, region3);
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(chimericAlignment);
        Assert.assertTrue(breakpoints.complication.getHomologyForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getInsertedSequenceForwardStrandRep(), "GAGATAGAGTC");
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for CIGAR extraction on tandem duplications
    // -----------------------------------------------------------------------------------------------
    @Test
    public void testExtractCigar() {

        final int contigTotalLength = 355;

        // forward strand
        final AlignmentRegion region1 = new AlignmentRegion(AlignmentRegion.DUMMY_ASM_ID, "1", new SimpleInterval("1", 1000001, 1000125),
                TextCigarCodec.decode("5H10S15M20D25M30D35M260S5H"),
                true, 60, 0, 16, 75);
        final AlignmentRegion region2 = new AlignmentRegion(AlignmentRegion.DUMMY_ASM_ID, "1", new SimpleInterval("1", 1000041, 1000145),
                TextCigarCodec.decode("5H185S45M30I55M20I5M10S5H"),
                true, 60, 0, 191, 340);

        final Cigar cigar1 = BreakpointComplications.extractCigarForTandup(region1, 1000125, 1000041);
        Assert.assertEquals(cigar1, TextCigarCodec.decode("20M30D35M"));
        final Cigar cigar2 = BreakpointComplications.extractCigarForTandup(region2, 1000125, 1000041);
        Assert.assertEquals(cigar2, TextCigarCodec.decode("45M30I40M"));

        // reverse strand
        final AlignmentRegion region3 = new AlignmentRegion(AlignmentRegion.DUMMY_ASM_ID, "1", region2.referenceInterval,
                CigarUtils.invertCigar(region2.cigarAlong5to3DirectionOfContig),
                false, 60, 0, contigTotalLength-region2.endInAssembledContig+1, contigTotalLength-region2.startInAssembledContig+1);
        final AlignmentRegion region4 = new AlignmentRegion(AlignmentRegion.DUMMY_ASM_ID, "1", region1.referenceInterval,
                CigarUtils.invertCigar(region1.cigarAlong5to3DirectionOfContig),
                false, 60, 0, contigTotalLength-region1.endInAssembledContig+1, contigTotalLength-region1.startInAssembledContig+1);

        final Cigar cigar3 = BreakpointComplications.extractCigarForTandup(region3, 1000125, 1000041);
        Assert.assertEquals(CigarUtils.invertCigar(cigar3), cigar2);
        final Cigar cigar4 = BreakpointComplications.extractCigarForTandup(region4, 1000125, 1000041);
        Assert.assertEquals(CigarUtils.invertCigar(cigar4), cigar1);
    }
}
package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.io.IOException;
import java.util.*;

/**
 * A helper struct for annotating complications that make the locations represented by its associated {@link NovelAdjacencyReferenceLocations}
 * a little ambiguous, so that downstream analysis could infer sv type with these complications.
 * To be updated as more types of complications can be processed and handled by
 * {@link #resolveComplications(ChimericAlignment)}.
 */
@DefaultSerializer(BreakpointComplications.Serializer.class)
final class BreakpointComplications {

    /**
     * In {@link BreakpointComplications#resolveComplications(ChimericAlignment)} where the naive attempt to resolve number of tandem repeats
     * on the reference and sample is done, we assume the lower number of repeats is no higher than this number.
     */
    private static final int MAX_LOWER_CN = 10;

    /**
     * '+' strand representations of micro-homology, inserted sequence and duplicated sequence on the reference.
     */
    final String homologyForwardStrandRep;
    final String insertedSequenceForwardStrandRep;
    final SimpleInterval dupSeqRepeatUnitRefSpan;
    final int dupSeqRepeatNumOnRef;
    final int dupSeqRepeatNumOnCtg;
    final List<String> cigarStringsForDupSeqOnCtg;
    final boolean dupAnnotFromOptimization;

    static final SimpleInterval DUPSEQ_REPEAT_UNIT_NA_VALUE
            = new SimpleInterval("NA", Integer.MAX_VALUE, Integer.MAX_VALUE);
    private static final BreakpointComplications UNHANDLED_CASE
            = new BreakpointComplications("", "", DUPSEQ_REPEAT_UNIT_NA_VALUE, 0, 0, new ArrayList<>(), false);


    /**
     * Given an ChimericAlignment representing two reference intervals rearranged as two intervals on the locally-assembled contig,
     * identify potential complications such as homology and duplication on the reference and/or on the contig.
     */
    @VisibleForTesting
    static BreakpointComplications resolveComplications(final ChimericAlignment chimericAlignment)
            throws IOException {

        final SimpleInterval leftReferenceInterval  = chimericAlignment.getCoordSortedReferenceIntervals()._1,
                             rightReferenceInterval = chimericAlignment.getCoordSortedReferenceIntervals()._2;

        final AlignmentRegion firstContigRegion  = chimericAlignment.regionWithLowerCoordOnContig,
                              secondContigRegion = chimericAlignment.regionWithHigherCoordOnContig;

        final byte[] contigSeq = chimericAlignment.contigSeq;

        // a segment with lower coordinate on the locally-assembled contig could map to a higher reference coordinate region
        // under two basic types of SV's: inversion (strand switch necessary) and translocation (no strand switch necessary)
        final boolean isNotSimpleTranslocation = ChimericAlignment.isNotSimpleTranslocation(chimericAlignment.regionWithLowerCoordOnContig, chimericAlignment.regionWithHigherCoordOnContig,
                chimericAlignment.strandSwitch, ChimericAlignment.involvesRefPositionSwitch(firstContigRegion, secondContigRegion));
        if (chimericAlignment.strandSwitch!= ChimericAlignment.StrandSwitch.NO_SWITCH) { // the case involves an inversion
            // TODO: 12/5/16 duplication detection to be done for inversion alleles
            @SuppressWarnings("unchecked")
            final BreakpointComplications result = new BreakpointComplications(getHomology(firstContigRegion, secondContigRegion, contigSeq),
                    getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq),
                    DUPSEQ_REPEAT_UNIT_NA_VALUE, 0, 0, Collections.EMPTY_LIST, false); // the "un-necessary" local var is for the warning suppression annotation
            return result;
        } else if (isNotSimpleTranslocation) {
            return getLocationComplicationForInsDel(chimericAlignment, leftReferenceInterval, rightReferenceInterval, contigSeq);
        } else { // TODO: 12/5/16 simple translocation, don't tackle yet
            return UNHANDLED_CASE;
        }
    }

    /**
     * Given an {@link ChimericAlignment} representing two reference intervals rearranged as two intervals on the locally-assembled contig,
     * identify potential complications such as homology and duplication on the reference and/or on the contig.
     */
    @SuppressWarnings("unchecked")
    BreakpointComplications(final ChimericAlignment chimericAlignment) {
        final SimpleInterval leftReferenceInterval  = chimericAlignment.getCoordSortedReferenceIntervals()._1,
                rightReferenceInterval = chimericAlignment.getCoordSortedReferenceIntervals()._2;

        final AlignmentRegion firstContigRegion  = chimericAlignment.regionWithLowerCoordOnContig,
                secondContigRegion = chimericAlignment.regionWithHigherCoordOnContig;

        final byte[] contigSeq = chimericAlignment.contigSeq;

        // a segment with lower coordinate on the locally-assembled contig could map to a higher reference coordinate region
        // under two basic types of SV's: inversion (strand switch necessary) and translocation (no strand switch necessary)
        final boolean isNotSimpleTranslocation = ChimericAlignment.isNotSimpleTranslocation(chimericAlignment.regionWithLowerCoordOnContig, chimericAlignment.regionWithHigherCoordOnContig,
                chimericAlignment.strandSwitch, ChimericAlignment.involvesRefPositionSwitch(firstContigRegion, secondContigRegion));

        if (chimericAlignment.strandSwitch!= ChimericAlignment.StrandSwitch.NO_SWITCH) { // the case involves an inversion
            // TODO: 12/5/16 duplication detection to be done for inversion alleles
            homologyForwardStrandRep = getHomology(firstContigRegion, secondContigRegion, contigSeq);
            insertedSequenceForwardStrandRep = getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq);
            dupSeqRepeatUnitRefSpan = DUPSEQ_REPEAT_UNIT_NA_VALUE;
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
            cigarStringsForDupSeqOnCtg = Collections.EMPTY_LIST;
            dupAnnotFromOptimization = false;
        } else if (isNotSimpleTranslocation) {
            final BreakpointComplications proxy = getLocationComplicationForInsDel(chimericAlignment, leftReferenceInterval, rightReferenceInterval, contigSeq);
            homologyForwardStrandRep = proxy.homologyForwardStrandRep;
            insertedSequenceForwardStrandRep = proxy.insertedSequenceForwardStrandRep;
            dupSeqRepeatUnitRefSpan = proxy.dupSeqRepeatUnitRefSpan;
            dupSeqRepeatNumOnRef = proxy.dupSeqRepeatNumOnRef;
            dupSeqRepeatNumOnCtg = proxy.dupSeqRepeatNumOnCtg;
            cigarStringsForDupSeqOnCtg = proxy.cigarStringsForDupSeqOnCtg;
            dupAnnotFromOptimization = proxy.dupAnnotFromOptimization;
        } else { // TODO: 12/5/16 simple translocation, don't tackle yet
            homologyForwardStrandRep = UNHANDLED_CASE.homologyForwardStrandRep;
            insertedSequenceForwardStrandRep = UNHANDLED_CASE.insertedSequenceForwardStrandRep;
            dupSeqRepeatUnitRefSpan = UNHANDLED_CASE.dupSeqRepeatUnitRefSpan;
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = UNHANDLED_CASE.dupSeqRepeatNumOnCtg;
            cigarStringsForDupSeqOnCtg = UNHANDLED_CASE.cigarStringsForDupSeqOnCtg;
            dupAnnotFromOptimization = UNHANDLED_CASE.dupAnnotFromOptimization;
        }
    }

    private static BreakpointComplications getLocationComplicationForInsDel(final ChimericAlignment chimericAlignment,
                                                                            final SimpleInterval leftReferenceInterval, final SimpleInterval rightReferenceInterval,
                                                                            final byte[] contigSeq) {

        final AlignmentRegion firstContigRegion  = chimericAlignment.regionWithLowerCoordOnContig,
                              secondContigRegion = chimericAlignment.regionWithHigherCoordOnContig;

        final int r1e = leftReferenceInterval.getEnd(),
                  r2b = rightReferenceInterval.getStart(),
                  c1e = firstContigRegion.endInAssembledContig,
                  c2b = secondContigRegion.startInAssembledContig;

        final int distBetweenAlignRegionsOnRef = r2b - r1e - 1, // distance-1 between the two regions on reference, denoted as d1 in the comments below
                  distBetweenAlignRegionsOnCtg = c2b - c1e - 1; // distance-1 between the two regions on contig, denoted as d2 in the comments below

        String homologyForwardStrandRepresentation="", insertedSeqForwardStrandRepresentation="";
        SimpleInterval dupSeqRepeatUnitRefSpan = DUPSEQ_REPEAT_UNIT_NA_VALUE;
        final List<String> dupSeqCigarStrings = new ArrayList<>();
        int dupSeqRepeatNumOnRef=0, dupSeqRepeatNumOnCtg=0;
        if (distBetweenAlignRegionsOnRef>0 && distBetweenAlignRegionsOnCtg==0) {        // Deletion: simple deletion, deleted sequence is [r1e+1, r2b-1] on the reference
            homologyForwardStrandRepresentation    = "";
            insertedSeqForwardStrandRepresentation = "";
            dupSeqRepeatUnitRefSpan                = DUPSEQ_REPEAT_UNIT_NA_VALUE;
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
        } else if (distBetweenAlignRegionsOnRef==0 && distBetweenAlignRegionsOnCtg>0) { // Insertion: simple insertion, inserted sequence is the sequence [c1e+1, c2b-1] on the contig
            homologyForwardStrandRepresentation    = "";
            insertedSeqForwardStrandRepresentation = getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq);
            dupSeqRepeatUnitRefSpan                = DUPSEQ_REPEAT_UNIT_NA_VALUE;
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
        } else if (distBetweenAlignRegionsOnRef==0 && distBetweenAlignRegionsOnCtg<0) { // Deletion: tandem repeat contraction where reference has two copies but one copy was deleted on the contig; duplicated sequence on reference are [r1e-|d2|+1, r1e] and [r2b, r2b+|d2|-1]
            homologyForwardStrandRepresentation    = getHomology(firstContigRegion, secondContigRegion, contigSeq);
            insertedSeqForwardStrandRepresentation = "";
            dupSeqRepeatUnitRefSpan                = new SimpleInterval(leftReferenceInterval.getContig(), r1e-(c1e-c2b), r1e);
            dupSeqRepeatNumOnRef      = 2;
            dupSeqRepeatNumOnCtg      = 1;
        } else if (distBetweenAlignRegionsOnRef<0 && distBetweenAlignRegionsOnCtg>=0) { // Insertion: tandem repeat expansion of reference bases [r1e-|d1|+1, r1e] to contig bases [c1e-|d1|+1, c1e] and [c2b, c2b+|d1|-1]
            // with optional inserted sequence [c1e+1, c2b-1] in between the two intervals on contig
            homologyForwardStrandRepresentation    = "";
            insertedSeqForwardStrandRepresentation = distBetweenAlignRegionsOnCtg==0 ? "" : getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq); // note this does not incorporate the duplicated reference sequence
            dupSeqRepeatUnitRefSpan                = new SimpleInterval(leftReferenceInterval.getContig(), r2b, r1e);
            dupSeqRepeatNumOnRef      = 1;
            dupSeqRepeatNumOnCtg      = 2;
            if (firstContigRegion.forwardStrand) {
                dupSeqCigarStrings.add( TextCigarCodec.encode(extractCigarForTandup(firstContigRegion, r1e, r2b)) );
                dupSeqCigarStrings.add( TextCigarCodec.encode(extractCigarForTandup(secondContigRegion, r1e, r2b)) );
            } else {
                dupSeqCigarStrings.add( TextCigarCodec.encode(CigarUtils.invertCigar(extractCigarForTandup(firstContigRegion, r1e, r2b))) );
                dupSeqCigarStrings.add( TextCigarCodec.encode(CigarUtils.invertCigar(extractCigarForTandup(secondContigRegion, r1e, r2b))) );
            }
        } else if (distBetweenAlignRegionsOnRef>0 && distBetweenAlignRegionsOnCtg>0) {  // Deletion: deletion with scar, i.e. large non-conserved substitution, reference bases [r1e+1, r2b-1] is substituted with contig bases [c1e+1, c2b-1]
            homologyForwardStrandRepresentation    = "";
            insertedSeqForwardStrandRepresentation = getInsertedSequence(firstContigRegion, secondContigRegion, contigSeq);
            dupSeqRepeatUnitRefSpan                = DUPSEQ_REPEAT_UNIT_NA_VALUE;
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
        } else if (distBetweenAlignRegionsOnRef>0 && distBetweenAlignRegionsOnCtg<0) {  // Deletion: a sequence of bases of length d1+HOM is deleted, and there's homology (which could be dup, but cannot tell): leftFlank+HOM+[r1e+1, r2b-1]+HOM+rightFlank -> leftFlank+HOM+rightFlank
            homologyForwardStrandRepresentation    = getHomology(firstContigRegion, secondContigRegion, contigSeq);
            insertedSeqForwardStrandRepresentation = "";
            dupSeqRepeatUnitRefSpan                = DUPSEQ_REPEAT_UNIT_NA_VALUE;
            dupSeqRepeatNumOnRef = dupSeqRepeatNumOnCtg = 0;
        } else if (distBetweenAlignRegionsOnRef<0 && distBetweenAlignRegionsOnCtg<0) {  // most complicated case, see below
            // Deletion:  duplication with repeat number N1 on reference, N2 on contig, such that N1 <= 2*N2 (and N2<N1);
            // Insertion: duplication with repeat number N1 on reference, N2 on contig, such that N2 <= 2*N1 (and N1<N2);
            // in both cases, the equal sign on the right can be taken only when there's pseudo-homology between starting bases of the duplicated sequence and starting bases of the right flanking region
            // the reference system with a shorter overlap (i.e. with less-negative distance between regions) has a higher repeat number

            final TandemRepeatComplication duplicationComplication = new TandemRepeatComplication(distBetweenAlignRegionsOnRef, distBetweenAlignRegionsOnCtg);

//            final byte[] dupSeq;
//            if (chimericAlignment.isForwardStrandRepresentation){
//                dupSeq = Arrays.copyOfRange(contigSeq, c2b-1, c2b+duplicationComplication.repeatedSeqLen-1);
//            } else {
//                dupSeq = Arrays.copyOfRange(contigSeq, c2b+duplicationComplication.pseudoHomologyLen -1, c2b+duplicationComplication.pseudoHomologyLen +duplicationComplication.repeatedSeqLen-1);
//                SequenceUtil.reverseComplement(dupSeq);
//            }

            final boolean isExpansion = distBetweenAlignRegionsOnRef<distBetweenAlignRegionsOnCtg;

            final int repeatUnitSpanStart = r1e - duplicationComplication.pseudoHomologyLen - duplicationComplication.repeatedSeqLen*duplicationComplication.lowerRepeatNumberEstimate + 1;
            final int repeatUnitSpanEnd   = repeatUnitSpanStart + duplicationComplication.repeatedSeqLen - 1;
            homologyForwardStrandRepresentation    = getHomology(firstContigRegion, secondContigRegion, contigSeq);
            insertedSeqForwardStrandRepresentation = "";
            dupSeqRepeatUnitRefSpan                = new SimpleInterval(firstContigRegion.referenceInterval.getContig(), repeatUnitSpanStart, repeatUnitSpanEnd);
            dupSeqRepeatNumOnRef      = isExpansion ? duplicationComplication.lowerRepeatNumberEstimate  : duplicationComplication.higherRepeatNumberEstimate;
            dupSeqRepeatNumOnCtg      = isExpansion ? duplicationComplication.higherRepeatNumberEstimate : duplicationComplication.lowerRepeatNumberEstimate;

            return new BreakpointComplications(homologyForwardStrandRepresentation, insertedSeqForwardStrandRepresentation, dupSeqRepeatUnitRefSpan, dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg, dupSeqCigarStrings, true);
        } else if (distBetweenAlignRegionsOnRef==0 && distBetweenAlignRegionsOnCtg==0) {// SNP & indel
            throw new GATKException("Detected badly parsed chimeric alignment for identifying SV breakpoints; no rearrangement found: " + chimericAlignment.toString());
        }

        if (insertedSeqForwardStrandRepresentation.isEmpty()){
            if (dupSeqRepeatNumOnCtg != dupSeqRepeatNumOnRef && dupSeqRepeatUnitRefSpan.equals(DUPSEQ_REPEAT_UNIT_NA_VALUE))
                throw new GATKException("An identified breakpoint pair seem to suggest insertion but the inserted sequence is empty: " + chimericAlignment.toString());
        }

        return new BreakpointComplications(homologyForwardStrandRepresentation, insertedSeqForwardStrandRepresentation, dupSeqRepeatUnitRefSpan, dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg, dupSeqCigarStrings, false);
    }

    /**
     * Given a {@link AlignmentRegion} signalling a tandem duplication, extract corresponding CIGAR between the suspected repeated
     * sequence as marked on reference between [{@code r2b}, {@code r1e}], and the sequence on the contig.
     */
    @VisibleForTesting
    static Cigar extractCigarForTandup(final AlignmentRegion contigRegion,
                                       final int r1e, final int r2b) {

        final List<CigarElement> elementList = contigRegion.cigarAlong5to3DirectionOfContig.getCigarElements();
        final List<CigarElement> result = new ArrayList<>(elementList.size());
        int refBasesConsumed = 0;
        boolean initiatedCollection = false;
        final int refStart = contigRegion.referenceInterval.getStart(),
                  refEnd = contigRegion.referenceInterval.getEnd();
        final boolean isForwardStrand = contigRegion.forwardStrand;
        if (isForwardStrand) {
            for(final CigarElement cigarElement : elementList) {
                final CigarOperator operator = cigarElement.getOperator();
                if ( !operator.isClipping() ) {
                    refBasesConsumed += operator.consumesReferenceBases() ? cigarElement.getLength() : 0;
                    if ( refStart+refBasesConsumed > r2b){ // entered suspected repeat region on ref
                        if (initiatedCollection){
                            if (refStart+refBasesConsumed <= r1e+1) {
                                result.add(cigarElement);
                            } else {
                                result.add(new CigarElement(cigarElement.getLength()-(refStart+refBasesConsumed-r1e)+1, CigarOperator.M));
                                break;
                            }
                        } else {
                            if (refStart+refBasesConsumed <= r1e+1) { // entering for the first time and the span doesn't overshoot
                                result.add(new CigarElement(refStart+refBasesConsumed-r2b, CigarOperator.M));
                            } else {                                // entering for the first time but the span overshoots pass the repeated region on ref
                                result.add(new CigarElement(cigarElement.getLength()-(refStart+refBasesConsumed-r1e)+1, CigarOperator.M));
                                break; // just one step overshoots
                            }
                            initiatedCollection = true;
                        }
                    }
                }
            }
        } else {
            for(final CigarElement cigarElement : elementList) {
                final CigarOperator operator = cigarElement.getOperator();
                if ( !operator.isClipping() ) {
                    refBasesConsumed += operator.consumesReferenceBases() ? cigarElement.getLength() : 0;
                    if ( refEnd-refBasesConsumed < r1e ) {
                        if (initiatedCollection){
                            if (refEnd-refBasesConsumed >= r2b-1) {
                                result.add(cigarElement);
                            } else {
                                result.add(new CigarElement(cigarElement.getLength()-(r2b-refEnd+refBasesConsumed)+1, CigarOperator.M));
                                break;
                            }
                        } else {
                            if (refEnd-refBasesConsumed >= r2b-1) {
                                result.add(new CigarElement(r1e - (refEnd-refBasesConsumed), CigarOperator.M));
                            } else {
                                result.add(new CigarElement(cigarElement.getLength()-(r2b-(refEnd-refBasesConsumed))+1, CigarOperator.M));
//                                result.add(new CigarElement(refBasesConsumed - (refEnd - r1e), CigarOperator.M));
                                break;
                            }
                            initiatedCollection = true;
                        }
                    }
                }
            }
        }
        return new Cigar(result);
    }

    /**
     * @return Micro-homology sequence using two alignments of the same contig: as indicated by their overlap on the contig itself.
     *          Empty if they don't overlap on the contig.
     */
    @VisibleForTesting
    static String getHomology(final AlignmentRegion current, final AlignmentRegion next, final byte[] contigSequence) {

        if (current.endInAssembledContig >= next.startInAssembledContig) {
            final byte[] homologyBytes = Arrays.copyOfRange(contigSequence, next.startInAssembledContig-1, current.endInAssembledContig);
            if (current.referenceInterval.getStart() > next.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(homologyBytes, 0, homologyBytes.length);
            }
            return new String(homologyBytes);
        } else {
            return "";
        }
    }

    /**
     * Note: not suitable for the most complicated case dealt with in {@link BreakpointComplications#resolveComplications(ChimericAlignment)}
     * @return Inserted sequence using two alignments of the same contig: as indicated by their separation on the the contig itself.
     */
    @VisibleForTesting
    static String getInsertedSequence(final AlignmentRegion current, final AlignmentRegion next, final byte[] contigSequence) {

        if (current.endInAssembledContig < next.startInAssembledContig-1) {
            final byte[] insertedSequenceBytes = Arrays.copyOfRange(contigSequence, current.endInAssembledContig, next.startInAssembledContig-1);
            if (current.referenceInterval.getStart() > next.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(insertedSequenceBytes, 0, insertedSequenceBytes.length);
            }
            return new String(insertedSequenceBytes);
        } else {
            return "";
        }
    }

    // TODO: 03/03/17 this complicated tandem duplication case is not exactly reproducible (e.g. '+' and '-' strand may give slightly different results by this treatment)
    private static final class TandemRepeatComplication {
        final int lowerRepeatNumberEstimate;
        final int higherRepeatNumberEstimate;
        final int repeatedSeqLen;
        final int pseudoHomologyLen;

        /**
         * This function, when given overlaps of two corresponding regions on reference and contig sequences,
         * attempts to find--naively and slowly--the repeat numbers on the reference and on the contig of tandem repeats,
         * as well as the pseudo-homology between the duplicated sequence and the right flanking region.
         *
         * An example might help:
         * an assembled contig that's actually a repeat expansion from 1 repeat to 2 repeats with pseudo-homology:
         * TGCCAGGTTACATGGCAAAGAGGGTAGATATGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
         * can be aligned to chr18:312XXX:
         * the 1st alignment chr18:312579-718 140M135S, which can be broken into the following part
         * 31:  TGCCAGGTTACATGGCAAAGAGGGTAGATAT
         * 109: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAA
         * 135: GAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
         * And the arithmetic works this way:
         * 31 + 109 = 140
         * 109 = 96 + 13
         * where 31 is the left flanking region before the repeated unit, which itself is 96 bases long (see below),
         * the number 13 is the length of the pseudo-homology between the starting bases of the repeated sequence and the right flanking region
         * a clearer picture emerges when we look at the 2nd alignment chr18:312610-757 127S148M, which can be broken into
         * 31: TGCCAGGTTACATGGCAAAGAGGGTAGATAT
         * 96: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA
         * 96: GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA
         * 13: GGGCAGCTGTGGA
         * 39: TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC
         * And the arithmetic works this way:
         * 31 + 96 = 127
         * 96 + 13 + 39 = 148
         */
        @VisibleForTesting
        TandemRepeatComplication(final int distBetweenAlignRegionsOnRef, final int distBetweenAlignRegionsOnCtg) {
            // the reference system with a shorter overlap (i.e. with less-negative distance between regions) has a higher repeat number
            final boolean isExpansion = distBetweenAlignRegionsOnRef<distBetweenAlignRegionsOnCtg;
            final int overlapOnLowerCNSequence, overlapOnHigherCNSequence;
            if (isExpansion) {
                overlapOnLowerCNSequence = Math.abs(distBetweenAlignRegionsOnRef);
                overlapOnHigherCNSequence = Math.abs(distBetweenAlignRegionsOnCtg);
            } else {     // d1 is lower absolute value -> reference has higher copy number of the duplication, i.e. Deletion
                overlapOnLowerCNSequence = Math.abs(distBetweenAlignRegionsOnCtg);
                overlapOnHigherCNSequence = Math.abs(distBetweenAlignRegionsOnRef);
            }

            int higherCnEst=0, lowerCnEst=0, unitLen=0, pseudoHomLen=0;
            double err = Double.MAX_VALUE;
            for(int cn2 = 1; cn2< MAX_LOWER_CN; ++cn2) {
                for(int cn1=cn2+1; cn1<=2*cn2; ++cn1) {
                    final int dupLenUpperBound = cn1==2*cn2 ? overlapOnLowerCNSequence : overlapOnHigherCNSequence;
                    for (int l=2; l<=dupLenUpperBound; ++l) {
                        for (int lambda=0; lambda<l; ++lambda) {
                            final int d1 = (2*cn2 - cn1)*l + lambda;
                            final int d2 = cn2*l + lambda;
                            final double newErr = Math.abs(overlapOnHigherCNSequence-d1) + Math.abs(overlapOnLowerCNSequence-d2);
                            if (newErr < err) {
                                err = newErr;
                                higherCnEst = cn1; lowerCnEst = cn2;
                                unitLen= l; pseudoHomLen = lambda;
                            }
                            if (err<1){
                                lowerRepeatNumberEstimate = lowerCnEst;
                                higherRepeatNumberEstimate = higherCnEst;
                                repeatedSeqLen = unitLen;
                                pseudoHomologyLen = pseudoHomLen;
                                return;
                            }
                        }
                    }
                }
            }

            lowerRepeatNumberEstimate = lowerCnEst;
            higherRepeatNumberEstimate = higherCnEst;
            repeatedSeqLen = unitLen;
            pseudoHomologyLen = pseudoHomLen;
        }
    }


    protected BreakpointComplications(final String homologyForwardStrandRep,
                                      final String insertedSequenceForwardStrandRep,
                                      final SimpleInterval dupSeqRepeatUnitRefSpan,
                                      final int dupSeqRepeatNumOnRef,
                                      final int dupSeqRepeatNumOnCtg,
                                      final List<String> cigarStringsForDupSeq,
                                      final boolean dupAnnotFromOptimization) {
        this.homologyForwardStrandRep = homologyForwardStrandRep;
        this.insertedSequenceForwardStrandRep = insertedSequenceForwardStrandRep;
        this.dupSeqRepeatUnitRefSpan = dupSeqRepeatUnitRefSpan;
        this.dupSeqRepeatNumOnRef = dupSeqRepeatNumOnRef;
        this.dupSeqRepeatNumOnCtg = dupSeqRepeatNumOnCtg;
        this.cigarStringsForDupSeqOnCtg = new ArrayList<>(cigarStringsForDupSeq);
        this.dupAnnotFromOptimization = dupAnnotFromOptimization;
    }

    protected BreakpointComplications(final Kryo kryo, final Input input) {
        homologyForwardStrandRep = input.readString();
        insertedSequenceForwardStrandRep = input.readString();
        dupSeqRepeatUnitRefSpan = new SimpleInterval(input.readString(), input.readInt(), input.readInt());
        dupSeqRepeatNumOnRef = input.readInt();
        dupSeqRepeatNumOnCtg = input.readInt();
        final int cigarCounts = input.readInt();
        cigarStringsForDupSeqOnCtg = new ArrayList<>(dupSeqRepeatNumOnCtg);
        for(int i=0; i<cigarCounts; ++i) {
            cigarStringsForDupSeqOnCtg.add(input.readString());
        }
        dupAnnotFromOptimization = input.readBoolean();
    }

    @Override
    public String toString() {
        return String.format("homology: %s\tinserted sequence: %s\ttandem duplication repeat unit ref span: %s\tref repeat num: %d\tctg repeat num: %d\tcigarStringsForDupSeqOnCtg: %s",
                homologyForwardStrandRep, insertedSequenceForwardStrandRep, dupSeqRepeatUnitRefSpan, dupSeqRepeatNumOnRef, dupSeqRepeatNumOnCtg, cigarStringsForDupSeqOnCtg);
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(homologyForwardStrandRep);
        output.writeString(insertedSequenceForwardStrandRep);
        output.writeString(dupSeqRepeatUnitRefSpan.getContig());
        output.writeInt(dupSeqRepeatUnitRefSpan.getStart());
        output.writeInt(dupSeqRepeatUnitRefSpan.getEnd());
        output.writeInt(dupSeqRepeatNumOnRef);
        output.writeInt(dupSeqRepeatNumOnCtg);
        output.writeInt(cigarStringsForDupSeqOnCtg.size());
        for(int i=0; i<cigarStringsForDupSeqOnCtg.size(); ++i) {
            output.writeString(cigarStringsForDupSeqOnCtg.get(i));
        }
        output.writeBoolean(dupAnnotFromOptimization);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        BreakpointComplications that = (BreakpointComplications) o;

        if (dupSeqRepeatNumOnRef != that.dupSeqRepeatNumOnRef) return false;
        if (dupSeqRepeatNumOnCtg != that.dupSeqRepeatNumOnCtg) return false;
        if (dupAnnotFromOptimization != that.dupAnnotFromOptimization) return false;
        if (!homologyForwardStrandRep.equals(that.homologyForwardStrandRep)) return false;
        if (!insertedSequenceForwardStrandRep.equals(that.insertedSequenceForwardStrandRep)) return false;
        if (!dupSeqRepeatUnitRefSpan.equals(that.dupSeqRepeatUnitRefSpan)) return false;
        return cigarStringsForDupSeqOnCtg.equals(that.cigarStringsForDupSeqOnCtg);
    }

    @Override
    public int hashCode() {
        int result = homologyForwardStrandRep.hashCode();
        result = 31 * result + insertedSequenceForwardStrandRep.hashCode();
        result = 31 * result + dupSeqRepeatUnitRefSpan.hashCode();
        result = 31 * result + dupSeqRepeatNumOnRef;
        result = 31 * result + dupSeqRepeatNumOnCtg;
        result = 31 * result + cigarStringsForDupSeqOnCtg.hashCode();
        result = 31 * result + (dupAnnotFromOptimization ? 1 : 0);
        return result;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<BreakpointComplications> {
        @Override
        public void write(final Kryo kryo, final Output output, final BreakpointComplications breakpointComplications) {
            breakpointComplications.serialize(kryo, output);
        }

        @Override
        public BreakpointComplications read(final Kryo kryo, final Input input, final Class<BreakpointComplications> klass ) {
            return new BreakpointComplications(kryo, input);
        }
    }
}

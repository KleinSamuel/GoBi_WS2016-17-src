package bamfiles;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

import augmentedTree.IntervalTree;
import genomeAnnotation.Exon;
import genomeAnnotation.Gene;
import genomeAnnotation.Transcript;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import util.Anders;
import util.GenRegVecUtil;
import util.Interval;

public class ReadPair {

	private LinkedList<Gene> spanningGenes = new LinkedList<>();
	private SAMRecord forward, reverse;
	private boolean isStarMapping;
	private LinkedList<Interval> blocksForward, blocksReverse;
	private LinkedList<Interval> intronsForward, intronsReverse;
	private LinkedList<Gene> matchedGenes = null; // stores transcript matches
													// or merged matches or
													// intronic matches
	private LinkedList<Transcript> matchedTranscripts;
	private int mismatchCount = 0, clippingSize = 0, splitCount = -2, geneDistance = Integer.MIN_VALUE, pcrIndex = 0;
	private boolean intergenic = false, antisense = false, intronic = false, merged = false, transcriptomic = false;
	private Interval genomicRegion;
	public String XX;
	private String calculatedXX = null;
	private boolean onNegativeStrand = false;

	public ReadPair(SAMRecord forward, SAMRecord reverse, boolean isStarMapping) {
		this.forward = forward;
		this.reverse = reverse;
		onNegativeStrand = forward.getReadNegativeStrandFlag();
		XX = (String) forward.getAttribute("XX");
		if (XX == null)
			XX = (String) reverse.getAttribute("XX");
		this.isStarMapping = isStarMapping;
		// blocksForward = GenRegVecUtil.parseAlignmentBlocks(forward);
		blocksForward = new LinkedList<>();
		for (AlignmentBlock ab : forward.getAlignmentBlocks()) {
			blocksForward.add(new Interval(ab.getReferenceStart() - 1, ab.getReferenceStart() + ab.getLength() - 2));
		}

		blocksForward = GenRegVecUtil.merge(blocksForward);
		intronsForward = GenRegVecUtil.getIntrons(blocksForward);
		blocksReverse = new LinkedList<>();
		for (AlignmentBlock ab : reverse.getAlignmentBlocks()) {
			blocksReverse.add(new Interval(ab.getReferenceStart() - 1, ab.getReferenceStart() + ab.getLength() - 2));
		}

		blocksReverse = GenRegVecUtil.merge(blocksReverse);
		// blocksReverse = GenRegVecUtil.parseAlignmentBlocks(reverse);
		intronsReverse = GenRegVecUtil.getIntrons(blocksReverse);
		genomicRegion = new Interval(Math.min(forward.getAlignmentStart(), reverse.getAlignmentStart()) - 1,
				Math.max(forward.getAlignmentEnd(), reverse.getAlignmentEnd()) - 1);
		calcSplitCount();
		if (splitCount == -1) {
			return;
		}

		if (splitCount > -1) {
			calcMissmatches();
			calcClipping();
			getMatchedGenes();
		}
	}

	public Interval getGenomicRegion() {
		return genomicRegion;
	}

	public boolean isOnNegativeStrand() {
		return onNegativeStrand;
	}

	public void calcMissmatches() {
		Integer missInForw = (Integer) forward.getAttribute("NM");
		Integer missInRev = (Integer) reverse.getAttribute("NM");
		if (missInForw == null || missInRev == null) {
			missInForw = (Integer) forward.getAttribute("nM");
			missInRev = (Integer) reverse.getAttribute("nM");
		}
		if (missInForw == null || missInRev == null) {
			missInForw = (Integer) forward.getAttribute("XM");
			missInRev = (Integer) reverse.getAttribute("XM");
		}
		if (isStarMapping) {
			mismatchCount = missInForw;
		} else {
			mismatchCount = missInForw + missInRev;
		}
	}

	public void calcClipping() {
		int clipp = Math.abs(forward.getAlignmentStart() - forward.getUnclippedStart())
				+ Math.abs(forward.getAlignmentEnd() - forward.getUnclippedEnd())
				+ Math.abs(reverse.getAlignmentStart() - reverse.getUnclippedStart())
				+ Math.abs(reverse.getAlignmentEnd() - reverse.getUnclippedEnd());

		this.clippingSize = clipp;
	}

	public int getSplitCount() {
		if (splitCount == -2) {
			calcSplitCount();
		}
		return splitCount;
	}

	/**
	 * returns -1 if splitInconsistent
	 */
	public void calcSplitCount() {
		if (checkIfSplitInconsistent()) {
			splitCount = -1;
		} else {
			splitCount = intronsForward.size() + intronsReverse.size();
			for (Interval i : intronsForward) {
				for (Interval j : intronsReverse) {
					if (i.overlaps(j)) {
						splitCount--;
					}

				}
			}
		}
	}

	public boolean checkIfSplitInconsistent() {
		for (Interval intronFw : intronsForward) {
			for (Interval i : blocksReverse) {
				if (i.overlaps(intronFw)) {
					return true;
				}
			}
		}
		for (Interval intronRv : intronsReverse) {
			for (Interval i : blocksForward) {
				if (i.overlaps(intronRv)) {
					return true;
				}
			}
		}
		return false;
	}

	public LinkedList<Gene> getSpanningGenes() {
		return spanningGenes;
	}

	public LinkedList<Gene> getMatchedGenes() {
		if (matchedGenes == null) {
			// get genes spanning forward read and genes spanning reverse reads
			// spanning means: gene.start <= start && gene.stop >= stop
			LinkedList<Gene> possibleGenes = Anders.ga.getChromosome(getReferenceName())
					.getSpecificStrandGenes(forward.getReadNegativeStrandFlag())
					.getIntervalsSpanning(genomicRegion.getStart(), genomicRegion.getStop(), new LinkedList<>());
			spanningGenes = possibleGenes;
			if (possibleGenes.size() == 0) {
				intergenic = Anders.ga.getChromosome(getReferenceName())
						.getSpecificStrandGenes(forward.getReadNegativeStrandFlag())
						.getIntervalsSpannedBy(genomicRegion.getStart(), genomicRegion.getStop(), new LinkedList<>())
						.isEmpty();
				if (!intergenic)
					geneDistance = 0;
				else
					getGeneDist();
				antisense = checkIfAntisense();
				matchedGenes = new LinkedList<>();
				return matchedGenes;
			}
			// for test:
			// matchedGenes = possibleGenes;
			// return matchedGenes;

			HashSet<Exon> nextExons = null, nextStrandExons = new HashSet<>();
			HashSet<Transcript> possibleTrsFW = new HashSet<>();// if
			// empty
			// -->
			// no
			// exon
			// spanning
			// regionBlock
			// -->
			// not
			// matching
			HashSet<Gene> forwardMatchingGenesMerged = new HashSet<>(); // matching
																		// at
																		// least
																		// merged
			HashSet<Gene> intronicGenes = new HashSet<>();
			HashSet<Transcript> possiblePerfectTrs = new HashSet<>();
			for (Gene g : possibleGenes) {
				HashSet<Transcript> trOverlaps = new HashSet<>();
				possibleTrsFW = new HashSet<>();
				for (Interval forwardInt : blocksForward) {
					nextExons = g.getAllExonsSorted().getIntervalsSpanning(forwardInt.getStart(), forwardInt.getStop(),
							new HashSet<>());
					if (!nextExons.isEmpty()) {
						nextStrandExons.addAll(nextExons);
						if (possibleTrsFW.isEmpty()) {
							for (Exon e : nextExons) {
								possibleTrsFW.addAll(e.getParentalTranscripts());
							}
						} else {
							for (Exon e : nextExons) {
								for (Transcript t : e.getParentalTranscripts()) {
									if (possibleTrsFW.contains(t)) {
										trOverlaps.add(t);
									}
								}
							}
							possibleTrsFW = new HashSet<>(trOverlaps);
							trOverlaps = new HashSet<>();
						}
					} else {
						nextStrandExons = new HashSet<>();
						trOverlaps = new HashSet<>();
						possibleTrsFW = new HashSet<>();
						intronicGenes.add(g);
						break;
					}
				}
				if (!nextStrandExons.isEmpty()) {
					forwardMatchingGenesMerged.add(g);
				}
				if (!possibleTrsFW.isEmpty()) {
					possiblePerfectTrs.addAll(possibleTrsFW);
				}
			}
			nextStrandExons = new HashSet<>();
			LinkedList<Gene> allMatched = new LinkedList<>();
			for (Gene g : forwardMatchingGenesMerged) {
				HashSet<Transcript> trOverlaps = new HashSet<>();
				for (Interval reverseInt : blocksReverse) {
					nextExons = g.getAllExonsSorted().getIntervalsSpanning(reverseInt.getStart(), reverseInt.getStop(),
							new HashSet<>());
					if (!nextExons.isEmpty()) {
						nextStrandExons.addAll(nextExons);
						if (!possiblePerfectTrs.isEmpty()) {
							for (Exon e : nextExons) {
								for (Transcript t : e.getParentalTranscripts()) {
									if (possiblePerfectTrs.contains(t)) {
										trOverlaps.add(t);
									}
								}
							}
							HashSet<Transcript> tmp = new HashSet<>();
							for (Transcript t : possiblePerfectTrs) {
								if (!t.getParentalGene().equals(g)) {
									tmp.add(t);
								}
							}
							possiblePerfectTrs = tmp;
							possiblePerfectTrs.addAll(trOverlaps);
							trOverlaps = new HashSet<>();
						}
					} else {
						nextStrandExons = new HashSet<>();
						intronicGenes.add(g);
						HashSet<Transcript> tmp = new HashSet<>();
						for (Transcript t : possiblePerfectTrs) {
							if (!t.getParentalGene().equals(g)) {
								tmp.add(t);
							}
						}
						possiblePerfectTrs = tmp;
						break;
					}
				}
				if (!nextStrandExons.isEmpty()) {
					allMatched.add(g);
				}
			}
			// all matched contains merged; possiblePerfectTrs contains possibly
			// perfect hit transcripts --> check if perfect hit
			// intronic genes contains intronic genes
			HashSet<Gene> matchedGenesAndTrs = new HashSet<>();
			matchedTranscripts = checkTranscripts(possiblePerfectTrs);
			for (Transcript t : matchedTranscripts) {
				matchedGenesAndTrs.add(t.getParentalGene());
			}
			matchedGenes = new LinkedList<>();
			matchedGenes.addAll(matchedGenesAndTrs);
			if (matchedGenes.isEmpty()) { // zwar keine transcript maps aber
											// merged
				for (Gene g : possibleGenes) {
					if (checkIfMerged(g)) {
						matchedGenes.add(g);
						merged = true;
					}
				}
				if (merged) {
					return matchedGenes;
				}
			} else {
				transcriptomic = true;
				return matchedGenes;
			}
			if (allMatched.isEmpty()) {
				// no transcript can map because merged is false --> intronic
				intronic = true;
				matchedGenes = new LinkedList<>();
				matchedGenes.addAll(intronicGenes);
				return matchedGenes;
			}
		}
		return matchedGenes;

	}

	public boolean checkIfAntisense() {
		LinkedList<Gene> possibleGenes = Anders.ga.getChromosome(getReferenceName())
				.getSpecificStrandGenes(!forward.getReadNegativeStrandFlag())
				.getIntervalsSpanning(genomicRegion.getStart(), genomicRegion.getStop(), new LinkedList<>());
		if (possibleGenes.isEmpty())
			return false;

		return true;

	}

	public boolean checkIfMerged(Gene g) {
		IntervalTree<Interval> mergedExons = g.getUnionTranscript();
		for (Interval forwardInt : blocksForward) {
			if (mergedExons.getIntervalsSpanning(forwardInt.getStart(), forwardInt.getStop(), new HashSet<>())
					.isEmpty())
				return false;
		}
		for (Interval reverseInt : blocksReverse) {
			if (mergedExons.getIntervalsSpanning(reverseInt.getStart(), reverseInt.getStop(), new HashSet<>())
					.isEmpty())
				return false;
		}
		return true;
	}

	public int getGeneDist() {
		if (geneDistance == Integer.MIN_VALUE) {
			LinkedList<Gene> neighbourLeft = Anders.ga.getChromosome(getReferenceName())
					.getSpecificStrandGenes(forward.getReadNegativeStrandFlag())
					.getIntervalsLeftNeighbor(genomicRegion.getStart(), genomicRegion.getStop(), new LinkedList<>()),
					neighbourRight = Anders.ga.getChromosome(getReferenceName())
							.getSpecificStrandGenes(forward.getReadNegativeStrandFlag()).getIntervalsRightNeighbor(
									genomicRegion.getStart(), genomicRegion.getStop(), new LinkedList<>());
			if (!neighbourLeft.isEmpty() || !neighbourRight.isEmpty()) {
				if (neighbourLeft.isEmpty()) {
					geneDistance = neighbourRight.getFirst().getStart() - genomicRegion.getStop() - 1;
				} else {
					if (neighbourRight.isEmpty()) {
						geneDistance = genomicRegion.getStart() - neighbourLeft.getFirst().getStop() - 1;
					} else {
						geneDistance = Math.min((genomicRegion.getStart() - neighbourLeft.getFirst().getStop() - 1),
								(neighbourRight.getFirst().getStart() - genomicRegion.getStop() - 1));
					}
				}
			}
			if (geneDistance < 0) {
				geneDistance = 0;
			}
		}
		return geneDistance;
	}

	private LinkedList<Transcript> checkTranscripts(HashSet<Transcript> possibleTrs) {
		LinkedList<Transcript> ret = new LinkedList<>();
		for (Transcript tr : possibleTrs) {
			if (checkTranscript(tr)) {
				ret.add(tr);
			}
		}
		return ret;
	}

	public boolean checkTranscript(Transcript tr) {
		boolean fw = false, rv = false;
		if (blocksForward.size() == 1) {
			fw = true;
		} else {
			for (Interval i : intronsForward) {
				if (tr.getIntrons().getIntervalsEqual(i.getStart(), i.getStop(), new LinkedList<>()).isEmpty())
					return false;
			}
			fw = true;
		}
		if (blocksReverse.size() == 1) {
			rv = true;
		} else {
			for (Interval i : intronsReverse) {
				if (tr.getIntrons().getIntervalsEqual(i.getStart(), i.getStop(), new LinkedList<>()).isEmpty())
					return false;
			}
			rv = true;
		}
		return fw && rv;
	}

	public String getAttributeXX() {
		if (calculatedXX == null) {
			calculatedXX = "mm:" + mismatchCount + "\tclipping:" + clippingSize + "\tgcount:" + matchedGenes.size()
					+ "\tnsplit:" + splitCount + "\t";
			if (matchedGenes.size() > 0) {
				Collections.sort(matchedGenes, new Comparator<Gene>() {

					@Override
					public int compare(Gene o1, Gene o2) {
						return o1.getId().compareTo(o2.getId());
					}
				});
				if (transcriptomic) {
					TreeMap<String, TreeSet<String>> sorted = new TreeMap<>();
					TreeSet<String> trs = null;
					for (Transcript t : matchedTranscripts) {
						trs = sorted.get(t.getParentalGene().getId());
						if (trs == null) {
							trs = new TreeSet<>();
							sorted.put(t.getParentalGene().getId(), trs);
						}
						trs.add(t.getId());
					}
					for (Entry<String, TreeSet<String>> e : sorted.entrySet()) {
						calculatedXX += e.getKey() + "," + Anders.ga.getGene(e.getKey()).getBiotype() + ":";
						for (String t : e.getValue()) {
							calculatedXX += t + ",";
						}
						calculatedXX = calculatedXX.substring(0, calculatedXX.length() - 1);
						calculatedXX += "\t";
					}
					calculatedXX = calculatedXX.substring(0, calculatedXX.length() - 1);
				} else {
					if (merged) {
						for (Gene g : matchedGenes) {
							calculatedXX += g.getId() + "," + g.getBiotype() + ":MERGED\t";
						}
						calculatedXX = calculatedXX.substring(0, calculatedXX.length() - 1);
					} else {
						if (intronic) {
							for (Gene g : matchedGenes) {
								calculatedXX += g.getId() + "," + g.getBiotype() + ":INTRON\t";
							}
							calculatedXX = calculatedXX.substring(0, calculatedXX.length() - 1);
						}
					}
				}
			} else {
				if (intergenic) {
					calculatedXX += "gdist:" + geneDistance + "\t";
				}
				calculatedXX += "antisense:" + antisense;
			}
			calculatedXX += "\tpcrindex: " + pcrIndex;
		}
		return calculatedXX;
	}

	public LinkedList<Interval> getForwardRegions() {
		return blocksForward;
	}

	public LinkedList<Interval> getReverseRegions() {
		return blocksReverse;
	}

	public SAMRecord getForward() {
		return forward;
	}

	public SAMRecord getReverse() {
		return reverse;
	}

	public String getReadName() {
		return forward.getReadName();
	}

	public String getReferenceName() {
		return forward.getReferenceName();
	}

	public LinkedList<Interval> getForwardIntrons() {
		return intronsForward;
	}

	public LinkedList<Interval> getReverseIntrons() {
		return intronsReverse;
	}

	public void setPCRindex(int index) {
		this.pcrIndex = index;
	}

	public LinkedList<Transcript> getMatchedTranscripts() {
		return matchedTranscripts;
	}

	public int getMismatchCount() {
		return mismatchCount;
	}

	public int getClippingSize() {
		return clippingSize;
	}

	public int getGeneDistance() {
		return geneDistance;
	}

	public int getPcrIndex() {
		return pcrIndex;
	}

	public boolean isIntergenic() {
		return intergenic;
	}

	public boolean isAntisense() {
		return antisense;
	}

	public boolean isIntronic() {
		return intronic;
	}

	public boolean isMerged() {
		return merged;
	}

	public boolean isTranscriptomic() {
		return transcriptomic;
	}

	public boolean compareXX() {
		return XX.equals(getAttributeXX());
	}

}

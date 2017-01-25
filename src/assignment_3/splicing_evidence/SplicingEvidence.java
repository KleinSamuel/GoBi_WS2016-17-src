package assignment_3.splicing_evidence;

import java.util.LinkedList;
import java.util.TreeSet;

import bamfiles.ReadPair;
import genomeAnnotation.Gene;
import genomeAnnotation.Transcript;
import util.Anders;
import util.Interval;

public class SplicingEvidence {

	private String chrId;
	private TreeSet<String> tr_ids;
	private TreeSet<String> geneId;
	private Interval intronRegion;
	private char strand = '+';
	private Boolean annotated = null;
	private int count = 0, ucount = 0, count_0 = 0, ucount_0 = 0, count_1 = 0, ucount_1 = 0, count_2 = 0, ucount_2 = 0;

	public SplicingEvidence(Interval intronRegion) {
		this.intronRegion = intronRegion;
		tr_ids = new TreeSet<>();
		geneId = new TreeSet<>();
		chrId = null;
	}

	public void addReadPair(ReadPair rp) {
		if (annotated == null) {
			checkIfAnnotated(rp);
		}
		updateEvidences(rp);
		if (rp.getPcrIndex() == 0) {
			ucount++;
			switch (rp.getMismatchCount()) {
			case 0:
				ucount_0++;
				ucount_1++;
				ucount_2++;
				break;
			case 1:
				ucount_1++;
				ucount_2++;
				break;
			case 2:
				ucount_2++;
				break;
			}
		}
		count++;
		switch (rp.getMismatchCount()) {
		case 0:
			count_0++;
			count_1++;
			count_2++;
			break;
		case 1:
			count_1++;
			count_2++;
			break;
		case 2:
			count_1++;
			break;
		}
	}

	public void checkIfAnnotated(ReadPair rp) {
		chrId = rp.getReferenceName();
		if (rp.isOnNegativeStrand())
			strand = '-';
		if (rp.isTranscriptomic()) {
			annotated = true;
		} else {
		LinkedList<Gene> possibleGenes = Anders.ga.getChromosome(chrId).getSpecificStrandGenes(strand == '-')
				.getIntervalsSpanning(getStart(), getStop(), new LinkedList<>());
		for (Gene g : possibleGenes) {
			for (Transcript t : g.getAllTranscriptsSorted()) {
				if (!t.getIntrons().getIntervalsEqual(getStart(), getStop(), new LinkedList<>()).isEmpty()) {
					annotated = true;
				}
			}
		}
		 }
		if (annotated == null) {
			annotated = false;
		}
		// wenn transcripts gematcht automatisch auch introns annotated, da
		// exon-bounds eingehalten sein müssen
		// if (rp.isTranscriptomic()) {
		// for (Gene g : rp.getMatchedGenes()) {
		// geneId.add(g.getId());
		// }
		// for (Transcript t : rp.getMatchedTranscripts()) {
		// tr_ids.add(t.getId());
		// }
		// return true;
		// }
		//
		// // wenn intergenic automatisch false, trivial
		// if (rp.isIntergenic()) {
		// return false;
		// }
		//
		// LinkedList<Gene> possibleGenes = null;
		// // wenn geneCount > 0 (merged oder intronic)
		// if (rp.isMerged() || rp.isIntronic()) {
		// possibleGenes = rp.getMatchedGenes();
		// } else { // wenn geneCount == 0, aber read intersected gen
		// possibleGenes =
		// Anders.ga.getChromosome(chrId).getSpecificStrandGenes(rp.isOnNegativeStrand())
		// .getIntervalsSpanning(getStart(), getStop(), new LinkedList<>());
		// }
		// for (Gene g : possibleGenes) {
		// for (Transcript t : g.getAllTranscriptsSorted()) {
		// if (t.getIntrons().contains(intronRegion)) {
		// geneId.add(g.getId());
		// tr_ids.add(t.getId());
		// }
		// }
		// }
		// if (geneId.isEmpty()) {
		// return false;
		// } else {
		// return true;
		// }
	}

	public void updateEvidences(ReadPair rp) {

		LinkedList<Gene> possibleGenes = Anders.ga.getChromosome(chrId).getSpecificStrandGenes(rp.isOnNegativeStrand())
				.getIntervalsSpanning(getStart(), getStop(), new LinkedList<>());
		for (Gene g : possibleGenes) {
			for (Transcript t : g.getAllTranscriptsSorted()) {
				if (!t.getIntrons().getIntervalsEqual(getStart(), getStop(), new LinkedList<>()).isEmpty()) {
					geneId.add(g.getId());
					tr_ids.add(t.getId());
				}
			}
		}

		// // wenn intergenic automatisch false, trivial
		// if (rp.isIntergenic()) {
		// return;
		// }
		//
		// // wenn transcripts gematcht automatisch auch introns annotated, da
		// // exon-bounds eingehalten sein müssen
		// if (rp.isTranscriptomic()) {
		// for (Gene g : rp.getMatchedGenes()) {
		// geneId.add(g.getId());
		// }
		// for (Transcript t : rp.getMatchedTranscripts()) {
		// tr_ids.add(t.getId());
		// }
		// return;
		// }
		//
		// LinkedList<Gene> possibleGenes = null;
		// // wenn geneCount > 0 (merged oder intronic)
		// if (rp.isMerged() || rp.isIntronic()) {
		// possibleGenes = rp.getMatchedGenes();
		// } else { // wenn geneCount == 0, aber read intersected gen
		// possibleGenes =
		// Anders.ga.getChromosome(chrId).getSpecificStrandGenes(rp.isOnNegativeStrand())
		// .getIntervalsSpanning(getStart(), getStop(), new LinkedList<>());
		// }
		// for (Gene g : possibleGenes) {
		// for (Transcript t : g.getAllTranscriptsSorted()) {
		// if (t.getIntrons().contains(intronRegion)) {
		// geneId.add(g.getId());
		// tr_ids.add(t.getId());
		// }
		// }
		// }

	}

	@Override
	public String toString() {
		String ret = chrId + "\t" + getStart() + "\t" + getStop() + "\t" + strand + "\t";
		if (geneId.size() > 0) {
			for (String gid : geneId) {
				ret += gid + "|";
			}
			ret = ret.substring(0, ret.length() - 1);
		}
		ret += "\t";
		if (tr_ids.size() > 0) {
			for (String trid : tr_ids) {
				ret += trid + "|";
			}
			ret = ret.substring(0, ret.length() - 1);
		}
		ret += "\t" + String.valueOf(annotated).toLowerCase() + "\t";
		ret += count + "\t" + ucount + "\t" + count_0 + "\t" + ucount_0 + "\t" + count_1 + "\t" + ucount_1 + "\t"
				+ count_2 + "\t" + ucount_2;
		return ret;
	}

	public String getChrId() {
		return chrId;
	}

	public TreeSet<String> getTr_ids() {
		return tr_ids;
	}

	public Boolean getAnnotated() {
		return annotated;
	}

	public int getCount() {
		return count;
	}

	public int getUcount() {
		return ucount;
	}

	public int getCount_0() {
		return count_0;
	}

	public int getUcount_0() {
		return ucount_0;
	}

	public int getCount_1() {
		return count_1;
	}

	public int getUcount_1() {
		return ucount_1;
	}

	public int getCount_2() {
		return count_2;
	}

	public int getUcount_2() {
		return ucount_2;
	}

	public int getStart() {
		return intronRegion.getStart();
	}

	public int getStop() {
		return intronRegion.getStop();
	}
}

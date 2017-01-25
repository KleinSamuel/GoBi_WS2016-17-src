package assignment_3;

import java.util.HashMap;

import bamfiles.ReadPair;
import genomeAnnotation.Transcript;

public class Task_3_Gene_Counts {

	private int NRPwithinGenReg, NURwithinGenReg = 0;
	private int NRPwithinGenRegWithoutOtherGenReg, NURwithinGenRegWithoutOtherGenReg = 0;
	private int NRPintronic, NURintronic = 0;
	private int NRPtranscriptomic, NURtranscriptomic = 0;
	private int NRPmergedTr, NURmergedTr = 0;
	private Integer NRPforHighestTranscript, NURforHighestTranscript = null;
	private HashMap<Transcript, Integer> NRPtrCounts, NURtrCounts;

	public Task_3_Gene_Counts() {
		this.NRPwithinGenReg = 0;
		this.NRPwithinGenRegWithoutOtherGenReg = 0;
		this.NRPintronic = 0;
		this.NRPtranscriptomic = 0;
		this.NRPmergedTr = 0;
		this.NRPforHighestTranscript = null;
		NRPtrCounts = new HashMap<>();
		NURtrCounts = new HashMap<>();
	}

	public void increaseNRPwithinGenReg(boolean NUR) {
		NRPwithinGenReg++;
		if (NUR) {
			NURwithinGenReg++;
		}
	}

	public void increaseNRPwithinGenRegWithout(boolean NUR) {
		NRPwithinGenRegWithoutOtherGenReg++;
		if (NUR) {
			NURwithinGenRegWithoutOtherGenReg++;
		}
	}

	public void addReadPair(ReadPair rp) {
		if (rp.isIntronic()) {
			NRPintronic++;
			if (rp.getPcrIndex() == 0) {
				NURintronic++;
			}
		} else {
			if (rp.isTranscriptomic()) {
				NRPtranscriptomic++;
				if (rp.getPcrIndex() == 0) {
					NURtranscriptomic++;
				}
			} else {
				if (rp.isMerged()) {
					NRPmergedTr++;
					if (rp.getPcrIndex() == 0) {
						NURmergedTr++;
					}
				}
			}
		}
	}

	public void addTranscript(Transcript t) {
		Integer i = NRPtrCounts.get(t);
		if (i == null) {
			NRPtrCounts.put(t, 1);
		} else {
			NRPtrCounts.put(t, i + 1);
		}
	}

	public void addNURTranscript(Transcript t) {
		Integer i = NURtrCounts.get(t);
		if (i == null) {
			NURtrCounts.put(t, 1);
		} else {
			NURtrCounts.put(t, i + 1);
		}
	}

	public int getNRPwithinGenReg() {
		return NRPwithinGenReg;
	}

	public int getNRPwithinGenRegWithoutOtherGenReg() {
		return NRPwithinGenRegWithoutOtherGenReg;
	}

	public int getNRPintronic() {
		return NRPintronic;
	}

	public int getNRPtranscriptomic() {
		return NRPtranscriptomic;
	}

	public int getNRPmergedTr() {
		return NRPmergedTr;
	}

	public int getNRPforHighestTranscript() {
		if (NRPforHighestTranscript == null) {
			NRPforHighestTranscript = 0;
			for (Integer i : NRPtrCounts.values()) {
				NRPforHighestTranscript = Math.max(NRPforHighestTranscript, i);
			}
			NRPtrCounts = null;
		}
		return NRPforHighestTranscript;
	}

	public int getNURforHighestTranscript() {
		if (NURforHighestTranscript == null) {
			NURforHighestTranscript = 0;
			for (Integer i : NURtrCounts.values()) {
				NURforHighestTranscript = Math.max(NURforHighestTranscript, i);
			}
			NURtrCounts = null;
		}
		return NURforHighestTranscript;
	}

	@Override
	public String toString() {
		return "" + NRPwithinGenReg + "\t" + NRPwithinGenRegWithoutOtherGenReg + "\t" + NRPintronic + "\t"
				+ NRPtranscriptomic + "\t" + NRPmergedTr + "\t" + getNRPforHighestTranscript() + "\t" + NURwithinGenReg
				+ "\t" + NURwithinGenRegWithoutOtherGenReg + "\t" + NURintronic + "\t" + NURtranscriptomic + "\t"
				+ NURmergedTr + "\t" + getNURforHighestTranscript();
	}

}

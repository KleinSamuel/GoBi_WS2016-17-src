package util;

import assignment_3.Task_3_Gene_Counts;
import genomeAnnotation.Gene;

public class FPKMUtils {

	public static double tr_fpkm(Gene g, Task_3_Gene_Counts geneCounts, double numberAllRPs) {
		double d = (Math.pow(10, 9) * geneCounts.getNRPtranscriptomic())
				/ (numberAllRPs * g.getMergedTranscript().getExonicLength());
		return d;
	}

	public static double intron_fpkm(Gene g, Task_3_Gene_Counts geneCounts, double numberAllRPs) {
		double d = (Math.pow(10, 9) * geneCounts.getNRPintronic())
				/ (numberAllRPs * g.getMergedTranscript().getExonicLength());
		return d;
	}

	public static double covP(Task_3_Gene_Counts geneCounts) {
		if (geneCounts.getNRPtranscriptomic() == 0) {
			return -1;
		}
		double d = 100d * geneCounts.getNRPforHighestTranscript() / geneCounts.getNRPtranscriptomic();
		return d;
	}

}

package bamfiles;

import java.io.BufferedWriter;
import java.util.Iterator;
import java.util.TreeMap;

public class PCRIndexStructure {

	// store genomic region vectors that implement comparable
	TreeMap<GenomicRegionVector, Integer> pcrIndicesPos;
	TreeMap<GenomicRegionVector, Integer> pcrIndicesNeg;
	String currentReferenceName = null;

	public PCRIndexStructure() {
		pcrIndicesPos = new TreeMap<>();
		pcrIndicesNeg = new TreeMap<>();
	}

	public int addReadPair(ReadPair rp, BufferedWriter bw) {
		if (currentReferenceName == null) {
			currentReferenceName = rp.getReferenceName();
		} else {
			if (!currentReferenceName.equals(rp.getReferenceName())) {
				currentReferenceName = rp.getReferenceName();
				pcrIndicesPos = new TreeMap<>();
				pcrIndicesNeg = new TreeMap<>();
			}
		}
		GenomicRegionVector grv = new GenomicRegionVector(rp.getForwardRegions(), rp.getReverseRegions());
		// delete some regions if they can't appear anymore
		if (pcrIndicesNeg.size() > 500000) {
			Iterator<GenomicRegionVector> iter = pcrIndicesNeg.keySet().iterator();

			while (iter.hasNext() && iter.next().getStart() < grv.getStart() / 2) {
				iter.remove();
			}
		}
		if (pcrIndicesPos.size() > 500000) {
			Iterator<GenomicRegionVector> iter = pcrIndicesPos.keySet().iterator();

			while (iter.hasNext() && iter.next().getStart() < grv.getStart() / 2) {
				iter.remove();
			}
		}

		try {
			// bw.write(grv.toString() + "\t" + rp.isOnNegativeStrand() + "\t" +
			// rp.getReferenceName() + "\t"
			// + rp.getReadName() + "\n");
		} catch (Exception e) {
			e.printStackTrace();
		}
		// System.out.println(grv);
		Integer pcrIndex = null;
		if (rp.isOnNegativeStrand()) {
			pcrIndex = pcrIndicesNeg.get(grv);
			if (pcrIndex == null) {
				pcrIndicesNeg.put(grv, 0);
				return 0;
			}
			pcrIndicesNeg.put(grv, pcrIndex + 1);
		} else {
			pcrIndex = pcrIndicesPos.get(grv);
			if (pcrIndex == null) {
				pcrIndicesPos.put(grv, 0);
				return 0;
			}
			pcrIndicesPos.put(grv, pcrIndex + 1);
		}
		return pcrIndex + 1;
	}

}

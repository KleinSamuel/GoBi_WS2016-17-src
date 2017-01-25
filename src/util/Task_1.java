package util;

import io.ConfigReader;
import io.HeadedFileReader;
import io.HeadedFileReader.LineObject;

public class Task_1 {

	public static void main(String[] args) {

		// read mappinginfo
		HeadedFileReader mappingInfo = new HeadedFileReader(ConfigReader.readConfig().get("mappinginfo_a3t2"), "\t");
		mappingInfo.readHeadedFile();
		String pathToEBVbams = ConfigReader.readConfig().get("EBVbams"), gtfPath = ConfigReader.readConfig().get("gtf");
		for (LineObject mapping : mappingInfo.getLineObjects()) {
			new Anders(pathToEBVbams + mapping.getValue("bamfile"), gtfPath, mapping).readBAMFile();
		}
	}
}

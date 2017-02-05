package io;

import java.io.File;
import java.io.IOException;

import debugStuff.DebugMessageFactory;

/**
 * Handle temporary files
 * 
 * @author Samuel Klein
 */
public class TemporaryFile {

	/**
	 * Create a temporary file which self destroys after the java VM terminates.
	 * 
	 * @param filepath
	 *            to directory for temporary file
	 * @return File temporary file
	 */
	public static File createTempFile(String filepath) {

		File tmp = null;

		try {
			tmp = File.createTempFile("tmp", ".txt", new File(filepath));
			// tmp.deleteOnExit();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return tmp;
	}

	/**
	 * Create a temporary file which self destroys after the java VM terminates.
	 * 
	 * @return File temporary file
	 */
	public static File createTempFile() {
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE, "Creating temp file at "+ConfigReader.readConfig().get("temp_directory"));
		return createTempFile(ConfigReader.readConfig().get("temp_directory"));
	}

}

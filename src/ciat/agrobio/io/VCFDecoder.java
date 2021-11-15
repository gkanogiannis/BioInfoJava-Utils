package ciat.agrobio.io;

import java.nio.ByteBuffer;

public class VCFDecoder implements VCFDecoderInterface<byte[]> {
	private static final byte LF = 10;
	//private boolean skip = false;

	public VCFDecoder() {
	}

	public byte[] decode(ByteBuffer buffer) {
		int lineStartPos = buffer.position();
		int limit = buffer.limit();
		while (buffer.hasRemaining()) {
			byte b = buffer.get();
			if (b == LF) { // reached line feed so parse line
				int lineEndPos = buffer.position();
				// set positions for one row duplication
				if (buffer.limit() < lineEndPos + 1) {
					buffer.position(lineStartPos).limit(lineEndPos);
				} else {
					buffer.position(lineStartPos).limit(lineEndPos + 1);
				}
				
				byte[] line = parseBytes(buffer, lineEndPos - lineStartPos - 1);
				
				if (line != null) {
					// reset main buffer
					buffer.position(lineEndPos);
					buffer.limit(limit);
					// set start after LF
					lineStartPos = lineEndPos;
				}
				if(line.length==0){
					continue;
				}
				
				/*
				if(line[0] == '>' || line[0] == '@'){
					skip = false;
				}
				else if(line[0] == '+'){
					skip = true;
				}
				*/
				
				/*
				if(skip){
					continue;
				}
				else{
					return line;
				}
				*/
				
				return line;
			}
		}
		buffer.position(lineStartPos);
		return null;
	}

	private byte[] parseBytes(ByteBuffer buffer, int length) {
		byte[] bytes = new byte[length];
		for (int i = 0; i < bytes.length; i++) {
			bytes[i] = buffer.get();
		}
		return bytes;
	}
}

interface VCFDecoderInterface<T> {

	public T decode(ByteBuffer buffer);

}
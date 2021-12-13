package runcycle;

public class State {
	private static NDDOConverter converter;

	public static NDDOConverter getConverter() {
		if (converter == null) throw new IllegalStateException("NDDOConverter not set!");

		return converter;
	}

	public static void setConverter(NDDOConverter converter) {
		if (converter == null) throw new NullPointerException("Converter cannot be null!");

		State.converter = converter;
	}
}

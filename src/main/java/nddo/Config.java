package nddo;

public class Config { // serializable
	public int rhf_numIt_max = 800;
	public int uhf_numIt_max = 1600;
	public double ediis_threshold = 1e-2, ediis_max_diff = 1e-5, diiserror_hard_limit = 1e-7;
	public double rhf_damp = 0.8, uhf_damp = 0.55;

	public double poplethiel_ideal = 1e-10;
	public int poplethiel_batch_size = 1;
	public int rhf_thiel_limit = 1000, uhf_thiel_limit = 1000;
	public int geom_opt_hessian_interval = 7;
}

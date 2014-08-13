source "tests_common.tcl"

require_feature "CUDA"

if { [catch { test_domain_decomposition_gpu } res] } {
	error_exit $res
} else {
	exit 0
}

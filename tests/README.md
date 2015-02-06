About
=====

Implementation of a small number of test cases.
Tests are implemented under the CuTest.h micro-framework and, as such, there are
various options for executions:
* Executions of the whole test set ( > tests)
* Execution of only some tests ( > tests [name_test1] [name_test2] ..)
* Executions skipping some test ( > tests --skip [name_test1] [name_test2] ..)

The suite contains the following tests:
* no_gaps_in_longer
* no_gaps_equal_length
* gaps_only_at_ends_in_shorter
* gaps_only_at_ends_in_longer
* no_gaps_smith_waterman
* free_start_gap
* free_end_gap
* free_gaps_at_ends
* free_gaps_at_ends_no_gaps
* free_gaps_at_ends_compare_no_gaps_gaps_only_at_ends
* no_mismatches_simple
* no_mismatches
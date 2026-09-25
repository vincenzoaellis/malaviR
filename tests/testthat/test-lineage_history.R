test_that("lineage_history_map reports unchanged, unknown, and follows a chain", {
  res <- lineage_history_map(c("SGS1", "NOT_A_LINEAGE_XYZ"))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("name", "current", "status", "trail"))
  expect_equal(res$current[1], "SGS1")
  expect_equal(res$status, c("unchanged", "unknown"))
  expect_true(is.na(res$current[2]))
  expect_equal(res$trail, c("", ""))
})

test_that("the walk itself follows renames, merges and removals", {
  ## The walk is exercised on a synthetic history through the same code path by
  ## temporarily substituting extract_table; the bundled release may have no
  ## history rows at all.
  fake_history <- data.frame(
    LINEAGE_NAME = c("A02", "B01", "T01", "S01"),
    EVENT        = c("renamed", "merged_into", "removed", "sequence_updated"),
    OLD_NAME     = c("A01", "A02", "T01", "S01"),
    NEW_NAME     = c("A02", "B01", "", "S01"),
    OLD_SEQUENCE = "", NEW_SEQUENCE = "", DATE = "", REASON = "", SOURCE = "", RELEASE = "",
    stringsAsFactors = FALSE)
  fake_gls <- data.frame(LINEAGE_NAME = c("B01", "S01", "U01"), stringsAsFactors = FALSE)
  with_mocked_bindings(
    extract_table = function(table, version = "latest") {
      if (table == "Lineage History") fake_history else fake_gls
    },
    {
      res <- lineage_history_map(c("A01", "A02", "T01", "S01", "U01", "Q99"))
      expect_equal(res$current, c("B01", "B01", NA, "S01", "U01", NA))
      ## A01 was renamed and then merged: its status is where it ended up
      expect_equal(res$status, c("merged", "merged", "removed", "sequence_updated",
                                 "unchanged", "unknown"))
      expect_equal(res$trail[1], "A01 renamed A02; A02 merged_into B01")
      expect_equal(res$trail[3], "T01 removed ")
    })
})

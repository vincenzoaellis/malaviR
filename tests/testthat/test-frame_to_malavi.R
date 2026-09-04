## frame_to_malavi places a primer-trimmed ASV into the 479 bp MalAvi frame.
## Built-in shortcuts (haem 2-479, leuc 2-477) and an arbitrary custom window
## must all pad the right number of N on each side and reject off-length input.

clean_seq <- function(n, base = "A") paste(rep(base, n), collapse = "")

test_that("haem shortcut pads one N on the left (478 -> 479, positions 2-479)", {
  framed <- frame_to_malavi(clean_seq(478), primer = "haem")
  expect_equal(nchar(framed), 479L)
  expect_equal(substr(framed, 1, 1), "N")
  expect_equal(substr(framed, 2, 2), "A")
  expect_equal(substr(framed, 479, 479), "A")   # no right pad
})

test_that("leuc shortcut pads one N left and two N right (476 -> 479, positions 2-477)", {
  framed <- frame_to_malavi(clean_seq(476), primer = "leuc")
  expect_equal(nchar(framed), 479L)
  expect_equal(substr(framed, 1, 1), "N")
  expect_equal(substr(framed, 478, 479), "NN")
})

test_that("a custom window (e.g. a ~456 bp set, positions 24-479) frames correctly", {
  framed <- frame_to_malavi(clean_seq(456), frame_start = 24, frame_end = 479)
  expect_equal(nchar(framed), 479L)
  expect_equal(substr(framed, 1, 23), clean_seq(23, "N"))   # 23 N on the left
  expect_equal(substr(framed, 24, 24), "A")
  expect_equal(substr(framed, 479, 479), "A")               # no right pad
})

test_that("custom window equal to the haem window matches the haem shortcut", {
  expect_identical(frame_to_malavi(clean_seq(478), frame_start = 2, frame_end = 479),
                   frame_to_malavi(clean_seq(478), primer = "haem"))
})

test_that("off-length ASVs are set to NA (default) with a warning, names preserved", {
  x <- c(good = clean_seq(478), bad = clean_seq(400))
  expect_warning(out <- frame_to_malavi(x, primer = "haem"))
  expect_equal(names(out), c("good", "bad"))
  expect_false(is.na(out[["good"]]))
  expect_true(is.na(out[["bad"]]))
})

test_that("on_off_length = 'keep' returns the original off-length sequence", {
  expect_warning(out <- frame_to_malavi(clean_seq(400), primer = "haem",
                                        on_off_length = "keep"))
  expect_equal(out, clean_seq(400))
})

test_that("a full-length (479 bp) sequence is treated as off-length, not reframed", {
  ## conservative behavior: an already-full-length sequence is not forced into
  ## the frame (it would shift the reading frame); default set_na returns NA
  expect_warning(out <- frame_to_malavi(clean_seq(479), primer = "haem"))
  expect_true(is.na(out))
})

test_that("on_off_length = 'error' stops on an off-length ASV", {
  expect_error(frame_to_malavi(clean_seq(400), primer = "haem",
                               on_off_length = "error"),
               "not the expected clean")
})

test_that("supplying both primer and an explicit window is an error", {
  expect_error(frame_to_malavi(clean_seq(478), primer = "haem",
                               frame_start = 2, frame_end = 479),
               "not both")
})

test_that("an out-of-range window is rejected", {
  expect_error(frame_to_malavi(clean_seq(10), frame_start = 0, frame_end = 479),
               "Invalid window")
  expect_error(frame_to_malavi(clean_seq(10), frame_start = 2, frame_end = 500),
               "Invalid window")
})

test_that("an unknown built-in primer name is rejected", {
  expect_error(frame_to_malavi(clean_seq(478), primer = "galen"),
               "must be one of")
})

test_that("an NA sequence is returned as NA instead of aborting the call", {
  ## nchar(NA) is NA, and an NA in a subscripted assignment is an error, so one
  ## empty cell in an ASV table used to kill the whole call with
  ## "NAs are not allowed in subscripted assignments".
  good <- paste(rep("A", 478), collapse = "")
  expect_warning(res <- frame_to_malavi(c(ok = good, missing = NA), primer = "haem"),
                 "are NA")
  expect_equal(nchar(res[["ok"]]), 479L)
  expect_true(is.na(res[["missing"]]))
  expect_equal(names(res), c("ok", "missing"))
})

test_that("an NA stays NA even when off-length sequences are kept", {
  good <- paste(rep("A", 478), collapse = "")
  short <- paste(rep("A", 400), collapse = "")
  res <- suppressWarnings(
    frame_to_malavi(c(ok = good, short = short, missing = NA),
                    primer = "haem", on_off_length = "keep"))
  expect_equal(res[["short"]], short)     # kept
  expect_true(is.na(res[["missing"]]))    # not kept, because it is not a sequence
})

test_that("interior whitespace is stripped, not treated as extra length", {
  ## The documentation says whitespace is stripped; trimws() only took the ends,
  ## so a 478 bp sequence with one internal space was reported off-length.
  good <- paste(rep("A", 478), collapse = "")
  spaced <- paste0(substr(good, 1, 200), " \n\t", substr(good, 201, 478))
  res <- frame_to_malavi(spaced, primer = "haem")
  expect_equal(nchar(res), 479L)
  expect_equal(res, frame_to_malavi(good, primer = "haem"))
})

test_that("pad_char must be a single character", {
  good <- paste(rep("A", 478), collapse = "")
  expect_error(frame_to_malavi(good, primer = "haem", pad_char = "NN"),
               "single character")
})

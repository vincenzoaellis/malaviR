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

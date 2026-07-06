## Tests for sister_taxa(): exercised on small, hand-checked ape trees so the
## expected output is transparent. Node numbering for these newick strings was
## confirmed directly from `tree$edge` (ape numbers tips 1..Ntip in order of
## appearance, then the root, then remaining internal nodes).

## ((A,B),(C,D));  ->  tips A=1 B=2 C=3 D=4, root=5, node 6=(A,B), node 7=(C,D)
two_cherries <- ape::read.tree(text = "((A,B),(C,D));")

## ((A,B),(C,(D,E)));  ->  node 8 = (C,(D,E)): a tip (C) beside a clade (D,E)
tip_plus_clade <- ape::read.tree(text = "((A,B),(C,(D,E)));")

## (((A,B),C),(D,E));  ->  the (((A,B),C)) node has a clade (A,B) beside a tip (C)
clade_plus_tip <- ape::read.tree(text = "(((A,B),C),(D,E));")

test_that("output has the documented columns", {
  res <- sister_taxa(two_cherries, node = 6)
  expect_s3_class(res, "data.frame")
  expect_identical(names(res), c("ancestral.node", "sister.clade", "taxa"))
})

test_that("a node subtending two tips returns both tips, one per clade", {
  res <- sister_taxa(two_cherries, node = 6)   # node 6 = (A, B)
  expect_equal(nrow(res), 2L)
  expect_equal(res$taxa, c("A", "B"))
  expect_equal(res$sister.clade, c(1, 2))
  expect_true(all(res$ancestral.node == 6))
})

test_that("a node subtending two clades returns both clades", {
  res <- sister_taxa(two_cherries, node = 5)   # root: (A,B) vs (C,D)
  expect_equal(nrow(res), 4L)
  ## clade 1 is the (A,B) side, clade 2 is the (C,D) side
  expect_setequal(res$taxa[res$sister.clade == 1], c("A", "B"))
  expect_setequal(res$taxa[res$sister.clade == 2], c("C", "D"))
})

test_that("a node subtending a tip and a clade returns both", {
  res <- sister_taxa(tip_plus_clade, node = 8) # node 8 = (C, (D,E))
  expect_equal(nrow(res), 3L)
  expect_setequal(res$taxa[res$sister.clade == 1], "C")
  expect_setequal(res$taxa[res$sister.clade == 2], c("D", "E"))
})

test_that("a node subtending a clade and a tip returns both (clade first)", {
  node <- ape::getMRCA(clade_plus_tip, c("A", "B", "C"))  # ((A,B),C)
  res <- sister_taxa(clade_plus_tip, node = node)
  expect_equal(nrow(res), 3L)
  expect_setequal(res$taxa[res$sister.clade == 1], c("A", "B"))
  expect_setequal(res$taxa[res$sister.clade == 2], "C")
})

test_that("multiple nodes are stacked into one data frame", {
  res <- sister_taxa(two_cherries, node = c(6, 7))
  expect_equal(nrow(res), 4L)                       # 2 tips per cherry
  expect_setequal(unique(res$ancestral.node), c(6, 7))
  expect_setequal(res$taxa[res$ancestral.node == 6], c("A", "B"))
  expect_setequal(res$taxa[res$ancestral.node == 7], c("C", "D"))
})

test_that("results are deterministic", {
  expect_identical(sister_taxa(tip_plus_clade, node = 8),
                   sister_taxa(tip_plus_clade, node = 8))
})

test_that("a tip number is rejected as not an internal node", {
  expect_error(sister_taxa(two_cherries, node = 1), "internal node")
})

test_that("an out-of-range node is rejected", {
  expect_error(sister_taxa(two_cherries, node = 99), "internal node")
})

test_that("non-tree input is rejected helpfully", {
  expect_error(sister_taxa(list(a = 1), node = 6), "class 'phylo'")
  expect_error(sister_taxa("not a tree", node = 6), "class 'phylo'")
})

test_that("a missing or non-numeric node is rejected", {
  expect_error(sister_taxa(two_cherries), "node")
  expect_error(sister_taxa(two_cherries, node = "6"), "node")
})

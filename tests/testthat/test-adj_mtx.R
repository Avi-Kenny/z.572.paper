
test_that("Adjacency matrix creator works (path)", {
  expect_equal(
    sum(adj_from_graph(generate_graph_path())),
    100 + (98*2 + 2*1)
  )
})

test_that("Adjacency matrix creator works (lattice)", {
  expect_equal(
    sum(adj_from_graph(generate_graph_lattice())),
    100 + (4*2 + 32*3 + 64*4)
  )
})

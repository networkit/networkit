networkit_add_module(matching
    BMatcher.cpp
    BMatching.cpp
    BSuitorMatcher.cpp
    LocalMaxMatcher.cpp
    Matcher.cpp
    Matching.cpp
    PathGrowingMatcher.cpp
    SuitorMatcher.cpp
    )

networkit_module_link_modules(matching
    auxiliary graph structures)

add_subdirectory(test)


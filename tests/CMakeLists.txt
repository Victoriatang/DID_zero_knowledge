add_executable(bbsac_bbs_test BBS_test.cpp)
target_link_libraries(bbsac_bbs_test PRIVATE bbsac_core)
add_test(NAME bbsac_bbs_test COMMAND $<TARGET_FILE:bbsac_bbs_test>)


add_executable(bbsac_range_proof_test range_proof_test.cpp)
target_link_libraries(bbsac_range_proof_test PRIVATE bbsac_core)
add_test(NAME bbsac_range_proof_test COMMAND $<TARGET_FILE:bbsac_range_proof_test>)

add_executable(bbsac_binarysearch_test binary_search.cpp)
target_link_libraries(bbsac_binarysearch_test PRIVATE bbsac_core)
add_test(NAME bbsac_binarysearch_test COMMAND $<TARGET_FILE:bbsac_binarysearch_test>)

add_executable(bbsac_issue_test issue_test.cpp)
target_link_libraries(bbsac_issue_test PRIVATE bbsac_core)
add_test(NAME bbsac_issue_test COMMAND $<TARGET_FILE:bbsac_issue_test>)

add_executable(bbsac_revocation_test revocation_test.cpp)
target_link_libraries(bbsac_revocation_test PRIVATE bbsac_core)
add_test(NAME bbsac_revocation_test COMMAND $<TARGET_FILE:bbsac_revocation_test>)


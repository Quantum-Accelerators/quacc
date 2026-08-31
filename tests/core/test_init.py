from __future__ import annotations

from logging import INFO, WARNING, getLogger


def test_emmet_task_doc_log_filter(caplog):
    logger = getLogger("emmet.core.tasks")

    with caplog.at_level(INFO, logger="emmet.core.tasks"):
        logger.info("Getting task doc in: /path/to/calculation")
        logger.info("Another useful Emmet message")
        logger.warning("Getting task doc in: warning context")

    assert "Getting task doc in: /path/to/calculation" not in caplog.text
    assert "Another useful Emmet message" in caplog.text
    assert "Getting task doc in: warning context" in caplog.text
    assert caplog.records[-1].levelno == WARNING

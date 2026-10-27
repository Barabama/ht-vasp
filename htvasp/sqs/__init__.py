"""SQS calculation building blocks with lazy package-level exports."""

from importlib import import_module


_EXPORTS = {
    "SqsCalculationMaker": ("htvasp.sqs.maker", "SqsCalculationMaker"),
    "GeneratedSqs": ("htvasp.sqs.generator", "GeneratedSqs"),
    "SqsGenerationError": ("htvasp.sqs.generator", "SqsGenerationError"),
    "SqsGenerator": ("htvasp.sqs.generator", "SqsGenerator"),
    "SqsgeneratorAdapter": ("htvasp.sqs.generator", "SqsgeneratorAdapter"),
    "SqsBatchManifest": ("htvasp.sqs.orchestrator", "SqsBatchManifest"),
    "SqsEndpointRecord": ("htvasp.sqs.orchestrator", "SqsEndpointRecord"),
    "SqsOrchestrationError": ("htvasp.sqs.orchestrator", "SqsOrchestrationError"),
    "SqsOrchestrator": ("htvasp.sqs.orchestrator", "SqsOrchestrator"),
    "SqsPointRecord": ("htvasp.sqs.orchestrator", "SqsPointRecord"),
    "SqsSystemRecord": ("htvasp.sqs.orchestrator", "SqsSystemRecord"),
    "endpoint_cache_key": ("htvasp.sqs.orchestrator", "endpoint_cache_key"),
    "SqsObservation": ("htvasp.sqs.task_doc", "SqsObservation"),
    "SqsStaticObservation": ("htvasp.sqs.task_doc", "SqsStaticObservation"),
    "make_static_observation": ("htvasp.sqs.jobs", "make_static_observation"),
    "normalize_static_task_document": ("htvasp.sqs.jobs", "normalize_static_task_document"),
    "sqs_observation": ("htvasp.sqs.jobs", "sqs_observation"),
}


def __getattr__(name: str):
    try:
        module_name, attribute = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
    value = getattr(import_module(module_name), attribute)
    globals()[name] = value
    return value

__all__ = [
    "SqsCalculationMaker",
    "GeneratedSqs",
    "SqsGenerationError",
    "SqsGenerator",
    "SqsgeneratorAdapter",
    "SqsBatchManifest",
    "SqsEndpointRecord",
    "SqsOrchestrationError",
    "SqsOrchestrator",
    "SqsSystemRecord",
    "SqsObservation",
    "SqsStaticObservation",
    "make_static_observation",
    "normalize_static_task_document",
    "sqs_observation",
    "endpoint_cache_key",
]

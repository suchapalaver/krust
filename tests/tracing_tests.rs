//! Tests for tracing instrumentation.
//!
//! These tests verify that tracing spans and events are emitted correctly
//! when the tracing feature is enabled.

#![cfg(feature = "tracing")]

use kmerust::run::count_kmers;
use kmerust::streaming::count_kmers_streaming;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use tracing::Level;
use tracing_subscriber::layer::SubscriberExt;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

/// A simple layer that counts events at INFO level or above.
struct EventCounter {
    count: Arc<AtomicUsize>,
}

impl<S: tracing::Subscriber> tracing_subscriber::Layer<S> for EventCounter {
    fn on_event(
        &self,
        event: &tracing::Event<'_>,
        _ctx: tracing_subscriber::layer::Context<'_, S>,
    ) {
        if event.metadata().level() <= &Level::INFO {
            self.count.fetch_add(1, Ordering::SeqCst);
        }
    }
}

#[test]
fn count_kmers_emits_tracing_events() {
    let event_count = Arc::new(AtomicUsize::new(0));
    let layer = EventCounter {
        count: Arc::clone(&event_count),
    };

    // Create a subscriber with our counting layer
    let subscriber = tracing_subscriber::registry().with(layer);

    // Run the test within the subscriber scope
    tracing::subscriber::with_default(subscriber, || {
        let path = fixture_path("simple.fa");
        let _counts = count_kmers(&path, 4).expect("should count k-mers");
    });

    // Should have emitted at least one event
    assert!(
        event_count.load(Ordering::SeqCst) > 0,
        "should emit tracing events"
    );
}

#[test]
fn streaming_count_emits_tracing_events() {
    let event_count = Arc::new(AtomicUsize::new(0));
    let layer = EventCounter {
        count: Arc::clone(&event_count),
    };

    let subscriber = tracing_subscriber::registry().with(layer);

    tracing::subscriber::with_default(subscriber, || {
        let path = fixture_path("simple.fa");
        let _counts = count_kmers_streaming(&path, 4).expect("should count k-mers");
    });

    assert!(
        event_count.load(Ordering::SeqCst) > 0,
        "should emit tracing events"
    );
}

#[test]
#[cfg(feature = "mmap")]
fn mmap_count_emits_tracing_events() {
    use kmerust::run::count_kmers_mmap;

    let event_count = Arc::new(AtomicUsize::new(0));
    let layer = EventCounter {
        count: Arc::clone(&event_count),
    };

    let subscriber = tracing_subscriber::registry().with(layer);

    tracing::subscriber::with_default(subscriber, || {
        let path = fixture_path("simple.fa");
        let _counts = count_kmers_mmap(&path, 4).expect("should count k-mers");
    });

    assert!(
        event_count.load(Ordering::SeqCst) > 0,
        "should emit tracing events"
    );
}
